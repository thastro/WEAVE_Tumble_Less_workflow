import os
import argparse


def setup_directories(base_dir="."):
    """
    Create the required WEAVE-TL XML directory and guide star catalogue structure
    under `base_dir`.

    Parameters
    ----------
    base_dir : str, optional
        Path where the structure will be created (default: current directory).
    """

    dirs = [
        "xmls",
        "xmls/plate_A/N",
        "xmls/plate_A/S",
        "xmls/plate_B/N",
        "xmls/plate_B/S",
        "xmls/configure_A/finalized_A",
        "xmls/configure_B/finalized_B",
        "resources/WG"

    ]

    for d in dirs:
        full_path = os.path.join(base_dir, d)
        os.makedirs(full_path, exist_ok=True)  # create if missing, no overwrite
        print(f"Ensured: {full_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create the WEAVE-TL XML directory structure."
    )
    parser.add_argument(
        "--base-dir", "-b",
        default=".",
        help="Base directory where the structure will be created (default: current directory)."
    )

    args = parser.parse_args()
    setup_directories(args.base_dir)
