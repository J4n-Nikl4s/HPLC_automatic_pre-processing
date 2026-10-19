import re
import pandas as pd
from pathlib import Path
from chromhandler import Handler, Molecule, Protein, to_enzymeml
from pathlib import Path
import numpy as np
import pyenzyme as pe

# unit used for concentration columns without a unit in the header
DEFAULT_CONC_UNIT = "umol / L"

# short unit notations and their chromhandler equivalent, other units are passed on unchanged
UNIT_ALIASES = {
    "M": "mol / L",
    "mM": "mmol / L",
    "uM": "umol / L",
    "µM": "umol / L",
    "μM": "umol / L",
    "nM": "nmol / L",
    "°C": "Celsius",
    "C": "Celsius",
}


def load_metadata(
    metadata_file: Path, sheet_name: str | int = 0
) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Load metadata from an Excel file and return a DataFrame with the metadata.

    Args:
        metadata_file (Path): Path to the Excel file containing the metadata.
        sheet_name (str | int): Name or index of the sheet. Defaults to the first sheet.

    Returns:
        tuple[pd.DataFrame, dict[str, str]]: DataFrame containing the metadata and
            the units from the column headers, e.g. {'TrisHCL': 'mmol / L'}.
    """
    metadata_df = pd.read_excel(metadata_file, sheet_name=sheet_name)


    # split column headers like 'TrisHCL [mM]' into name ('TrisHCL') and unit ('mmol / L')
    column_units = {}
    column_names = []
    for column in metadata_df.columns:
        match = re.fullmatch(r"(.+?)\s*\[(.+)\]", str(column).strip())
        if match:
            name, unit = match.group(1), match.group(2).strip()
            column_units[name] = UNIT_ALIASES.get(unit, unit)
        else:
            name = str(column).strip()
        column_names.append(name)
    metadata_df.columns = column_names

    # empty dilution factor cells mean the sample was not diluted
    if "dilution_factor" in metadata_df.columns:
        metadata_df["dilution_factor"] = metadata_df["dilution_factor"].fillna(1)

    metadata_df = metadata_df.dropna(axis='columns')

    # set 'experiment_name' as index
    metadata_df.set_index("experiment_name", inplace=True)

    print(f"Units from header: {column_units}")
    return metadata_df, column_units

def load_moleculecules(molecule_dirs: list[Path], retention_times_updated: dict[str, float]) -> list:
    """
    Load calibration molecules from JSON files in the specified directories. 
    Update their retention times based on the provided dictionary.
    Args:
        molecule_dirs (list[Path]): List of directories containing JSON files for calibration molecules.
        retention_times_updated (dict[str, float]): Dictionary mapping molecule IDs to their updated retention times.
    returns:
        list: List of loaded Molecule objects with updated retention times.
    """
    # loop over all the json files in the calibration molecules folders
    # and create molecule objects
    molecules_by_id = {}
    for molecule_dir in molecule_dirs:
        if not molecule_dir.is_dir():
            raise FileNotFoundError(f"Calibration molecule folder '{molecule_dir}' does not exist.")

        for molecule_path in molecule_dir.glob("*.json"):
            molecule = Molecule.read_json(molecule_path)
            molecules_by_id[molecule.id] = molecule

    # update retention times
    for molecule_id, retention_time in retention_times_updated.items():
        if molecule_id not in molecules_by_id:
            print(f"Warning: no molecule '{molecule_id}' loaded, retention time not updated.")
            continue
        molecules_by_id[molecule_id].retention_time = retention_time

    molecules = list(molecules_by_id.values())

    print(f"Loaded {len(molecules)} molecules:")
    for molecule in molecules:
        print(f"  {molecule.id}: retention time {molecule.retention_time} min")

    return molecules

def load_proteins(protein_dir: Path) -> list:
    """
    Load proteins from JSON files in the specified directory.

    Args:
        protein_dir (Path): Directory containing JSON files for proteins.

    Returns:
        list: List of loaded Protein objects.
    """
    proteins = []
    for protein_path in protein_dir:
        protein = Protein.read_json(protein_path)
        proteins.append(protein)

    print(f"Loaded {len(proteins)} with names {[p.name for p in proteins]}")
    return proteins
    
    
def find_id(name: str) -> str | None:
    """Returns the first dotted number in a name, e.g. '47.1.1_0.1 µM_EbPPK' -> '47.1.1'."""
    match = re.search(r"\d+(?:\.\d+)*", name)
    return match.group() if match else None


def get_metadata(sample_name: str, metadata_df: pd.DataFrame) -> pd.Series:
    """Returns the metadata row of a sample.

    Tries, in this order:
    1. a row with the exact sample name
    2. a row with the sample name without trailing replicate number, e.g. 'Assay 6.1' -> 'Assay 6'
    3. a row with the same ID as the sample, e.g. '47.1.1_...' -> 'Assay 47.1.1'
    4. a row with the experiment ID of the sample, e.g. '47.1.1_...' -> 'Assay 47.1'
    """
    if sample_name in metadata_df.index:
        return metadata_df.loc[sample_name]

    # strip the replicate number ('.1', ' 1', '_1', '-1') from the end of the name
    experiment_name = re.sub(r"[\s._-]*\d+$", "", sample_name)
    if experiment_name in metadata_df.index:
        return metadata_df.loc[experiment_name]

    # match by ID, the last number of the sample ID is the replicate number
    sample_id = find_id(sample_name)
    if sample_id is not None:
        row_ids = {row: find_id(str(row)) for row in metadata_df.index}
        candidate_ids = [sample_id]
        if "." in sample_id:
            candidate_ids.append(sample_id.rsplit(".", 1)[0])

        for candidate_id in candidate_ids:
            rows = [row for row, row_id in row_ids.items() if row_id == candidate_id]
            if len(rows) > 1:
                raise KeyError(
                    f"Sample '{sample_name}' matches several metadata rows: {rows}."
                )
            if rows:
                return metadata_df.loc[rows[0]]

    raise KeyError(
        f"No metadata for sample '{sample_name}'. Add a row named "
        f"'{sample_name}' or '{experiment_name}' or with ID '{sample_id}' to the metadata sheet."
    )

def analyze_time_samples(sample_dir: Path, molecules: list[Molecule], proteins: list[Protein], metadata_df: pd.DataFrame, column_units: dict[str, str]) -> list[Handler]:
    """
    Analyze all samples in the specified directory using the provided molecules and proteins.

    Args:
        sample_dir (Path): Directory containing sample data.
        molecules (list[Molecule]): List of Molecule objects to be used in the analysis.
        proteins (list[Protein]): List of Protein objects to be used in the analysis.
        metadata_df (pd.DataFrame): DataFrame containing metadata for each sample.
        column_units (dict[str, str]): Dictionary mapping column names to their units.

    Returns:
        list[Handler]: List of Handler objects containing the analysis results for each sample.
    """
    # loop over all the time samples and create the ChromAnalyzer objects
    handlers = []
    for time_samples in sample_dir.iterdir():
        path = str(time_samples)
        metadata = get_metadata(time_samples.name, metadata_df)

        analyzer = Handler.read_agilent(
            path=path,
            ph=metadata["pH"],
            temperature=metadata["temperature"],
            temperature_unit=column_units.get("temperature", "Celsius"),
            mode="timecourse",
        )

        # add molecule data to the analyzer
        for molecule in molecules:
            if molecule.id not in metadata_df.columns:
                continue
            init_conc = metadata[molecule.id]

            # set constant molecules (concentration does not change over time)
            if molecule.id in ["MgCl2", "TrisHCL"]:
                molecule.constant = True

            # if init conc not nan, add to analyzer
            if not np.isnan(init_conc):
                analyzer.add_molecule(
                    molecule,
                    init_conc=init_conc,
                    conc_unit=column_units.get(molecule.id, DEFAULT_CONC_UNIT),
                )

        # add protein data to the analyzer
        for protein in proteins:
            if protein.name not in metadata_df.columns:
                continue
            init_conc = metadata[protein.name]

            # if init conc not nan, add to analyzer
            if not np.isnan(init_conc):
                analyzer.add_protein(
                    protein,
                    init_conc=init_conc,
                    conc_unit=column_units.get(protein.name, DEFAULT_CONC_UNIT),
                )

        analyzer.assign_all_peaks()

        # set dilution factor (calculated concentrations are multiplied by it)
        if "dilution_factor" in metadata_df.columns:
            analyzer.set_dilution_factor(float(metadata["dilution_factor"]))

        handlers.append(analyzer)

    return handlers

def ensure_small_molecule(doc: pe.EnzymeMLDocument, id: str, name: str, **kwargs) -> None:
    """Adds a small molecule to the EnzymeML document if no species with this ID exists yet.

    Use it for species of a reaction that are not measured, e.g. PolyP.
    Further keyword arguments (e.g. canonical_smiles) are passed to add_to_small_molecules.
    """
    if id in _species_ids(doc):
        return
    vessel_id = doc.vessels[0].id if doc.vessels else None
    doc.add_to_small_molecules(id=id, name=name, vessel_id=vessel_id, **kwargs)


def _species_ids(doc: pe.EnzymeMLDocument) -> dict[str, str]:
    """Returns a mapping from species ID and name to species ID for all species of the document."""
    ids = {}
    for species in [*doc.small_molecules, *doc.proteins, *doc.complexes]:
        if species.name:
            ids[species.name] = species.id
        ids[species.id] = species.id
    return ids


def add_reaction(
    doc: pe.EnzymeMLDocument,
    id: str,
    name: str,
    reactants: dict[str, float],
    products: dict[str, float],
    modifiers: dict[str, str] | None = None,
    reversible: bool = False,
) -> pe.Reaction:
    """Adds a reaction to the EnzymeML document.

    Species can be given by ID or by name, e.g. the protein 'EbPPK' is stored with its ID 'P34367'.

    Args:
        doc: EnzymeML document to add the reaction to.
        id: ID of the reaction, e.g. 'r1'.
        name: Name of the reaction.
        reactants: species and their stoichiometry, e.g. {'AMP': 1, 'PolyP': 1}.
        products: species and their stoichiometry, e.g. {'ADP': 1}.
        modifiers: species and their role, e.g. {'EbPPK': 'BIOCATALYST'}. Roles are
            BIOCATALYST, CATALYST, ACTIVATOR, INHIBITOR, BUFFER, SOLVENT and ADDITIVE.
        reversible: whether the reaction is reversible.

    Returns:
        pe.Reaction: the added reaction.
    """
    modifiers = modifiers or {}
    ids = _species_ids(doc)

    # check all species before changing the document
    missing = [s for s in [*reactants, *products, *modifiers] if s not in ids]
    if missing:
        raise KeyError(
            f"Species {missing} of reaction '{id}' are not in the EnzymeML document. "
            f"Available: {sorted(set(ids.values()))}. Add missing species with ensure_small_molecule()."
        )
    if id in [reaction.id for reaction in doc.reactions]:
        raise ValueError(f"A reaction with ID '{id}' already exists in the EnzymeML document.")

    reaction = doc.add_to_reactions(id=id, name=name, reversible=reversible)
    for species, stoichiometry in reactants.items():
        reaction.add_to_reactants(species_id=ids[species], stoichiometry=stoichiometry)
    for species, stoichiometry in products.items():
        reaction.add_to_products(species_id=ids[species], stoichiometry=stoichiometry)
    for species, role in modifiers.items():
        role = role if isinstance(role, pe.ModifierRole) else pe.ModifierRole[str(role).upper()]
        reaction.add_to_modifiers(species_id=ids[species], role=role)

    return reaction
