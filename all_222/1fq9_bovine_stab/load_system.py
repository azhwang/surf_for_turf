import MDAnalysis as mda
import yaml


def load_config(filename="config.yaml"):

    with open(filename) as f:
        return yaml.safe_load(f)


def load_universe(config):

    u = mda.Universe(
        config["system"]["psf"],
        config["system"]["trajectory"]
    )

    return u


if __name__ == "__main__":

    cfg = load_config()

    u = load_universe(cfg)

    print("Atoms:", len(u.atoms))
    print("Frames:", len(u.trajectory))

    print("\nChains detected:")
    for chain in set(u.atoms.chainIDs):
        atoms = u.select_atoms(f"chainID {chain}")
        print(
            chain,
            atoms.residues.resnames[:5],
            "atoms:",
            len(atoms)
        )
