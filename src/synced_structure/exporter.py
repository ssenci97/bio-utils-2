"""Curated exporters for SyncedStructure: PDB and mmCIF structure files, TSV tables for chains, residues, moieties, alignments, contacts, and interfaces, single or grid figures, and standalone 3D view pages."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

import io
from collections.abc import Callable, Sequence
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from Bio.PDB import PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from tqdm import tqdm

from ._types import ChainPair, ChainRefs, Moieties, PathLike
from .base import StructureState
from .constants import (
    CONTACT_SEQUENCE_COLUMNS, CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES, DEFAULT_CONTACT_TYPE, EXPORT_SEP,
    EXPORT_SUFFIX, MAPPING_COLUMNS, MOIETY_PROTEIN, PDB_MAX_ATOMS, PDB_MAX_RESNUM, PDB_MIN_RESNUM,
    STRUCTURE_FORMATS,
)
from .panels import build_panel_grid, rasterize_figure
from .plot_style import EXPORT_DPI, FIGURE_FORMATS, PANEL_FORMATS, PDF_RC
from .viz3d import VIZ_PAGE_TEMPLATE

if TYPE_CHECKING:
    from Bio.PDB.Atom import Atom, DisorderedAtom
    from matplotlib.figure import Figure

    from .chains import ChainRecord

DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DI", "DU"}


def _atom_variants(atom: Atom | DisorderedAtom) -> list[Atom]:
    """Return every alternate location of an atom, or the atom itself.

    Args:
        atom: Biopython atom, disordered or not.

    Returns:
        The atoms to write.
    """
    if atom.is_disordered():
        return list(cast("DisorderedAtom", atom).disordered_get_list())
    return [cast("Atom", atom)]


def _cif(value: object) -> str:
    """Return an mmCIF value: '?' for missing, the text otherwise.

    Args:
        value: Value, possibly None or empty.

    Returns:
        The value as text.
    """
    text = "" if value is None else str(value).strip()
    return text or "?"


def _poly_type(record: ChainRecord) -> str:
    """Return the mmCIF entity_poly type of a polymer chain.

    Args:
        record: Chain record with polymer residues.

    Returns:
        'polypeptide(L)', 'polydeoxyribonucleotide', or 'polyribonucleotide'.
    """
    if record.moieties.count(MOIETY_PROTEIN) * 2 >= len(record.moieties):
        return "polypeptide(L)"
    names = [residue.get_resname().strip().upper() for residue in record.residues]
    dna = sum(name in DNA_RESIDUES for name in names)
    return "polydeoxyribonucleotide" if dna * 2 >= len(names) else "polyribonucleotide"

###################################### EXPORT HELPERS


class ExportMixin(StructureState):
    """Write tables and figures to the output directory."""

    def _resolve_outdir(self, outdir: PathLike | None = None) -> Path:
        """Resolve the output directory, creating it if missing.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The directory.
        """
        path = Path(outdir or self.config["save_outdir"]).expanduser()
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _write_table(self, frame: pd.DataFrame, name: str, outdir: PathLike | None = None) -> Path:
        """Write a table to '{structure_id}_{name}.tsv' in the output directory.

        Args:
            frame: Table to write.
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path.
        """
        path = self._resolve_outdir(outdir) / f"{self.structure_id}_{name}{EXPORT_SUFFIX}"
        frame.to_csv(path, sep=EXPORT_SEP, index=False)
        return path

    def _save_figure(self, fig: Figure, stem: str, outdir: PathLike | None, formats: Sequence[str],
                     dpi: int = EXPORT_DPI, rasterize_pdf: bool = True) -> list[Path]:
        """Save a figure at its exact size in each format, then close it.

        Args:
            fig: Figure to save.
            stem: File name without extension.
            outdir: Output directory, or None for the configured save_outdir.
            formats: File formats, such as 'pdf' or 'png'.
            dpi: Resolution of raster output.
            rasterize_pdf: Whether to rasterize collections and images in PDFs.

        Returns:
            The written paths.
        """
        base = self._resolve_outdir(outdir) / stem
        paths = []
        for fmt in formats:
            path = base.with_name(f"{base.name}.{fmt}")
            if fmt == "pdf" and rasterize_pdf:
                rasterize_figure(fig)
            # PDF_RC holds valid rc keys; matplotlib types them as a Literal union that a plain dict can't match
            with mpl.rc_context(cast("dict[Any, Any]", PDF_RC)):
                fig.savefig(path, dpi=dpi)
            paths.append(path)
        plt.close(fig)
        return paths

###################################### CHAIN AND RESIDUE TABLES

    def export_chain_table(self, outdir: PathLike | None = None) -> Path:
        """Export chain_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_chains.tsv'.
        """
        return self._write_table(self.chain_table(), "chains", outdir)

    def export_chain_residue_table(self, chain_ids: ChainRefs | None = None, outdir: PathLike | None = None) -> Path:
        """Export the residue-to-reference mapping of the selected chains.

        Args:
            chain_ids: Chains to export, or None for all.
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_residue_mapping.tsv'.
        """
        selected = self.chain_ids if chain_ids is None else [self._chain_key(c) for c in chain_ids]
        frames = [self.chain_residue_table(c) for c in tqdm(selected, desc="mapping", unit="chain", leave=False)]
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=MAPPING_COLUMNS)
        return self._write_table(frame, "residue_mapping", outdir)

    def export_residue_table(self, outdir: PathLike | None = None) -> Path:
        """Export residue_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_residues.tsv'.
        """
        return self._write_table(self.residue_table(), "residues", outdir)

    def export_moiety_table(self, outdir: PathLike | None = None) -> Path:
        """Export moiety_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_moieties.tsv'.
        """
        return self._write_table(self.moiety_table(), "moieties", outdir)

###################################### ALIGNMENT TABLES

    def export_alignment_table(self, outdir: PathLike | None = None) -> Path:
        """Export alignment_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_alignment_stats.tsv'.
        """
        return self._write_table(self.alignment_table(), "alignment_stats", outdir)

    def export_alignment_column_table(self, outdir: PathLike | None = None) -> Path:
        """Export alignment_column_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_alignment_positions.tsv'.
        """
        return self._write_table(self.alignment_column_table(), "alignment_positions", outdir)

###################################### CONTACT TABLES

    def export_contact_strings(self, outdir: PathLike | None = None, cutoff: float = CONTACTS_CUTOFF_NM,
                               contact_type: str = DEFAULT_CONTACT_TYPE, reference: bool = False,
                               scheme: str = "ca", moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> Path:
        """Export each chain's sequence and contact string.

        Args:
            outdir: Output directory, or None for the configured save_outdir.
            cutoff: Distance cutoff in nm.
            contact_type: Contact type.
            reference: Whether to use reference sequences (reference_contact_string), skipping chains without one.
            scheme: Contact scheme ('ca' or 'sidechain').
            moieties: Biomolecule types included in the output.

        Returns:
            The written path, '{structure_id}_contact_strings.tsv', or '{structure_id}_contact_strings_reference.tsv'
            with reference=True.
        """
        contacts = self.contact_table(cutoff=cutoff, contact_type=contact_type, scheme=scheme, moieties=moieties)
        options: dict[str, Any] = {"cutoff": cutoff, "contact_type": contact_type, "contacts": contacts,
                                   "scheme": scheme, "moieties": moieties}
        rows = []
        for chain_key, data in self.chains.items():
            if reference:
                if not data.reference_sequence:
                    continue
                sequence = data.reference_sequence
                text = self.reference_contact_string(chain_key, **options)
            else:
                sequence = self._moiety_sequence(chain_key, moieties)
                text = self.contact_string(chain_key, **options)
            rows.append({
                "structure_id": self.structure_id,
                "chain_key": chain_key,
                "chain_id": data.chain_id,
                "auth_chain_id": data.auth_chain_id,
                "sequence": sequence,
                "contacts": text,
                "contact_type": contact_type,
                "mdtraj_scheme": scheme,
                "cutoff_nm": cutoff,
            })
        name = "contact_strings_reference" if reference else "contact_strings"
        return self._write_table(pd.DataFrame(rows, columns=CONTACT_SEQUENCE_COLUMNS), name, outdir)

    def export_contact_stat_tables(self, outdir: PathLike | None = None, cutoff: float = CONTACTS_CUTOFF_NM,
                                   moieties: Moieties = DEFAULT_CONTACT_MOIETIES, scheme: str = "ca") -> list[Path]:
        """Export every table of contact_stat_tables.

        Args:
            outdir: Output directory, or None for the configured save_outdir.
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.

        Returns:
            The written paths, '{structure_id}_contact_stats_{name}.tsv' for pairs, profile, separation, and
            moiety_pairs.
        """
        tables = self.contact_stat_tables(cutoff=cutoff, moieties=moieties, scheme=scheme)
        return [self._write_table(frame, f"contact_stats_{name}", outdir) for name, frame in tables.items()]

###################################### INTERFACE TABLES

    def export_interface_tables(self, chain_pairs: Sequence[ChainPair] | None = None,
                                outdir: PathLike | None = None) -> list[Path]:
        """Export complex-level and residue-level interface metrics.

        Args:
            chain_pairs: Chain pairs to export, or None for all pairs.
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written paths, '{structure_id}_interface_complex.tsv' and '{structure_id}_interface_residues.tsv'.
        """
        tables = self.interface_tables(chain_pairs)
        return [self._write_table(tables["complex"], "interface_complex", outdir),
                self._write_table(tables["residues"], "interface_residues", outdir)]

###################################### FIGURES

    def export_plot(self, plot_func: Callable[[], tuple[Figure, Any]], name: str, outdir: PathLike | None = None,
                    formats: Sequence[str] = FIGURE_FORMATS, dpi: int = EXPORT_DPI, raster: bool = True) -> list[Path]:
        """Save the figure returned by plot_func() as '{structure_id}_{name}' in each format.

        Args:
            plot_func: Callable returning (figure, axes).
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.
            formats: File formats.
            dpi: Resolution of raster output.
            raster: Whether to rasterize collections and images in PDFs.

        Returns:
            The written paths.
        """
        fig, _ = plot_func()
        return self._save_figure(fig, f"{self.structure_id}_{name}", outdir, formats, dpi, raster)

    def export_plot_grid(self, plot_funcs: Sequence[Callable[..., Any]], name: str, outdir: PathLike | None = None,
                         formats: Sequence[str] = PANEL_FORMATS, dpi: int = EXPORT_DPI, raster: bool = True,
                         **grid_kwargs: Any) -> list[Path]:
        """Save a panel grid of plot functions as '{structure_id}_{name}' in each format.

        For several grids in one multi-page PDF, use export_panels_pdf.

        Args:
            plot_funcs: Plot functions, one per panel, as taken by build_panel_grid.
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.
            formats: File formats.
            dpi: Resolution of raster output.
            raster: Whether to rasterize collections and images in PDFs.
            **grid_kwargs: Passed to build_panel_grid.

        Returns:
            The written paths.
        """
        fig = build_panel_grid(plot_funcs, **grid_kwargs)
        return self._save_figure(fig, f"{self.structure_id}_{name}", outdir, formats, dpi, raster)

###################################### STRUCTURE FILES

    def _export_keys(self, chain_ids: ChainRefs | None) -> list[str]:
        """Return the chain keys to write, in model order.

        Args:
            chain_ids: Chains to write, or None for all.

        Returns:
            Chain keys.

        Raises:
            ValueError: If no chain is selected.
        """
        wanted = set(self.chain_ids if chain_ids is None else (self._chain_key(c) for c in chain_ids))
        keys = [key for key in self.chain_ids if key in wanted]
        if not keys:
            raise ValueError("No chains selected to write.")
        return keys

    def _entity_ids(self, keys: Sequence[str]) -> dict[str, str]:
        """Return the entity ID written for each chain: its own, unless chains sharing it differ in reference sequence.

        Args:
            keys: Chain keys.

        Returns:
            {chain key: entity ID}, one entity per (entity ID, reference sequence) pair.
        """
        used = {str(self.chains[key].entity_id) for key in keys if self.chains[key].entity_id}
        counter = 1
        owners: dict[tuple[str | None, str], str] = {}
        taken: set[str] = set()
        entities: dict[str, str] = {}
        for key in keys:
            record = self.chains[key]
            pair = (record.entity_id, record.reference_sequence)
            if pair not in owners:
                entity = record.entity_id
                if entity is None or entity in taken:
                    while str(counter) in used or str(counter) in taken:
                        counter += 1
                    entity = str(counter)
                owners[pair] = entity
                taken.add(entity)
            entities[key] = owners[pair]
        return entities

    def mmcif_string(self, chain_ids: ChainRefs | None = None) -> str:
        """Return the structure as mmCIF text that from_path() reads back with the same chains.

        Atoms keep their label and author chain IDs, and polymer residues their reference position as label_seq_id.
        Entities carry the reference sequence (entity_poly), description, and accession (struct_ref), so reference
        sequences, entity IDs, and accessions survive a round trip; chains sharing an entity ID but holding different
        reference sequences get separate entities. Metadata and refseq IDs set by the user are not written.

        Args:
            chain_ids: Chains to write, or None for all.

        Returns:
            The mmCIF text.
        """
        return self._mmcif_text(chain_ids)

    def _mmcif_text(self, chain_ids: ChainRefs | None = None, key_as_author: bool = False) -> str:
        """Return mmCIF text, see mmcif_string.

        Args:
            chain_ids: Chains to write, or None for all.
            key_as_author: Whether to write chain keys as author chain IDs, so viewers that name chains by author ID
                (3Dmol.js) show the keys; files written this way do not round-trip.

        Returns:
            The mmCIF text.
        """
        keys = self._export_keys(chain_ids)
        entities = self._entity_ids(keys)
        descriptions = self._entity_descriptions()
        atoms: dict[str, list[str]] = {name: [] for name in (
            "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
            "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
            "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num")}
        serial = 0
        for key in keys:
            record = self.chains[key]
            positions = {id(residue): record.sequence_index_to_ref_index.get(index)
                         for index, residue in enumerate(record.residues)}
            polymer = set(positions)
            for residue in record.chain:
                hetero = id(residue) not in polymer
                position = positions.get(id(residue))
                for atom in residue.get_atoms():
                    for variant in _atom_variants(atom):
                        serial += 1
                        name, resname = variant.get_id().strip(), residue.get_resname().strip()
                        for column, value in (
                                ("group_PDB", "HETATM" if hetero else "ATOM"), ("id", serial),
                                ("type_symbol", variant.element or "?"), ("label_atom_id", name),
                                ("label_alt_id", (variant.get_altloc() or "").strip() or "."),
                                ("label_comp_id", resname), ("label_asym_id", record.chain_id),
                                ("label_entity_id", entities[key]),
                                ("label_seq_id", "." if hetero or position is None else position + 1),
                                ("pdbx_PDB_ins_code", residue.id[2].strip() or "?"),
                                ("Cartn_x", f"{variant.coord[0]:.3f}"), ("Cartn_y", f"{variant.coord[1]:.3f}"),
                                ("Cartn_z", f"{variant.coord[2]:.3f}"), ("occupancy", f"{variant.get_occupancy():.2f}"),
                                ("B_iso_or_equiv", f"{variant.get_bfactor():.2f}"), ("auth_seq_id", residue.id[1]),
                                ("auth_comp_id", resname),
                                ("auth_asym_id", key if key_as_author else record.auth_chain_id),
                                ("auth_atom_id", name), ("pdbx_PDB_model_num", 1)):
                            atoms[column].append(str(value))
        data: dict[str, Any] = {"data_": self.structure_id}
        order = list(dict.fromkeys(entities[key] for key in keys))
        first = {entity: next(key for key in keys if entities[key] == entity) for entity in order}
        data["_entity.id"] = order
        data["_entity.type"] = ["polymer" if self.chains[first[e]].residues else "non-polymer" for e in order]
        data["_entity.pdbx_description"] = [_cif(descriptions.get(self.chains[first[e]].entity_id or "")) for e in order]
        polymers = [e for e in order if self.chains[first[e]].reference_sequence]
        if polymers:
            data["_entity_poly.entity_id"] = polymers
            data["_entity_poly.type"] = [_poly_type(self.chains[first[e]]) for e in polymers]
            for name in ("pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"):
                data[f"_entity_poly.{name}"] = [self.chains[first[e]].reference_sequence for e in polymers]
        asym = list(dict.fromkeys((self.chains[key].chain_id, entities[key]) for key in keys))
        data["_struct_asym.id"] = [label for label, _ in asym]
        data["_struct_asym.entity_id"] = [entity for _, entity in asym]
        refs = [(e, self.get_accession(first[e])) for e in order]
        refs = [(e, accession.strip()) for e, raw in refs if raw for accession in raw.split(",")]
        if refs:
            data["_struct_ref.id"] = [str(index) for index in range(1, len(refs) + 1)]
            data["_struct_ref.entity_id"] = [e for e, _ in refs]
            data["_struct_ref.pdbx_db_accession"] = [accession for _, accession in refs]
        data.update({f"_atom_site.{column}": values for column, values in atoms.items()})
        writer = MMCIFIO()
        writer.set_dict(data)
        output = io.StringIO()
        writer.save(output)
        return output.getvalue()

    def pdb_string(self, chain_ids: ChainRefs | None = None) -> str:
        """Return the structure as PDB text, with author chain IDs.

        PDB keeps coordinates, residue IDs, and author chain IDs only: reload an mmCIF file to keep reference
        sequences, entities, and label IDs.

        Args:
            chain_ids: Chains to write, or None for all.

        Returns:
            The PDB text.

        Raises:
            ValueError: If an author chain ID is not one character, two chains would merge, or atoms or residue
                numbers exceed the PDB format; write mmCIF instead.
        """
        keys = self._export_keys(chain_ids)
        structure = Structure(self.structure_id)
        model = Model(0)
        structure.add(model)
        atoms = 0
        for key in keys:
            record = self.chains[key]
            chain_id = record.auth_chain_id
            if len(chain_id) != 1:
                raise ValueError(f"Chain {key} has author ID {chain_id!r}; PDB needs one character. "
                                 "Rename it with rename_chains() or write mmCIF.")
            if chain_id not in model:
                model.add(Chain(chain_id))
            target = model[chain_id]
            for residue in record.chain:
                number = residue.id[1]
                if not PDB_MIN_RESNUM <= number <= PDB_MAX_RESNUM:
                    raise ValueError(f"Residue number {number} of chain {key} does not fit the PDB format; write mmCIF.")
                if residue.id in target:
                    raise ValueError(f"Chains sharing author ID {chain_id} both have residue {residue.id}; "
                                     "write mmCIF instead.")
                target.add(residue.copy())
                atoms += sum(len(_atom_variants(atom)) for atom in residue.get_atoms())
        if atoms > PDB_MAX_ATOMS:
            raise ValueError(f"{atoms} atoms exceed the PDB limit of {PDB_MAX_ATOMS}; write mmCIF.")
        writer = PDBIO()
        writer.set_structure(structure)
        output = io.StringIO()
        writer.save(output)
        return output.getvalue()

    def export_structure(self, file_format: str = "cif", chain_ids: ChainRefs | None = None,
                         outdir: PathLike | None = None, name: str | None = None) -> Path:
        """Write the structure, or some chains, to a PDB or mmCIF file.

        mmCIF (see mmcif_string) reloads with from_path() with the same chain keys, reference sequences, entities,
        and accessions; PDB (see pdb_string) keeps coordinates and author chain IDs only.

        Args:
            file_format: 'cif' (also 'mmcif') or 'pdb'.
            chain_ids: Chains to write, or None for all.
            outdir: Output directory, or None for the configured save_outdir.
            name: File name without extension, or None for the structure ID, followed by the chain keys when
                chain_ids is given.

        Returns:
            The written path, such as '1AAY.cif' or '1AAY_A-B.pdb'.

        Raises:
            ValueError: If file_format is not supported.
        """
        fmt = str(file_format).lower().lstrip(".")
        if fmt not in STRUCTURE_FORMATS:
            raise ValueError(f"file_format must be one of {sorted(STRUCTURE_FORMATS)}.")
        if name is None:
            name = str(self.structure_id)
            if chain_ids is not None:
                name += "_" + "-".join(self._export_keys(chain_ids))
        text = self.pdb_string(chain_ids) if STRUCTURE_FORMATS[fmt] == "pdb" else self.mmcif_string(chain_ids)
        path = self._resolve_outdir(outdir) / f"{name}.{STRUCTURE_FORMATS[fmt]}"
        path.write_text(text)
        return path

###################################### 3D VIEWS

    def export_structure_3d(self, name: str = "structure_3d", outdir: PathLike | None = None, title: str | None = None,
                            **view_kwargs: Any) -> Path:
        """Save show_structure_3d as a standalone HTML page '{structure_id}_{name}.html'.

        Args:
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.
            title: Optional title shown above the view and used as the page title.
            **view_kwargs: Passed to show_structure_3d.

        Returns:
            The written path.
        """
        view = self.show_structure_3d(title=title, **view_kwargs)
        page = VIZ_PAGE_TEMPLATE.replace("__TITLE__", escape(str(title or self.structure_id)))
        path = self._resolve_outdir(outdir) / f"{self.structure_id}_{name}.html"
        path.write_text(page.replace("__BODY__", view.data), encoding="utf-8")
        return path
