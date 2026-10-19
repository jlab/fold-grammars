import argparse
import json
import logging
import re
import subprocess
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests

logger: logging.Logger = logging.getLogger("RNAmotiFold_sequence_updater")

non_listed_conversions: dict[str, str] = {"MAD": "A"}
conversion_json_path: Path = (
    Path(__file__).resolve().parent.joinpath("nucleotide_conversion.json")
)


class MotifSequence:
    def __init__(
        self,
        motif_name: str = "",
        abbreviation: str = "",
        loop_type: str = "",
        file: str = "",
    ):
        self.motif_name: str = motif_name
        self.abbreviation: str = abbreviation
        self.loop_type: str = loop_type
        self.file: str = file

    @classmethod
    def from_dict(cls, motif_json_entry: dict[str, str]):
        """Alternative constructor to create MotifSequence object from RNA3D json files"""
        return cls(
            motif_name=motif_json_entry["motif_name"],
            abbreviation=motif_json_entry["abbreviation"],
            loop_type=motif_json_entry["loop_type"],
            file=motif_json_entry["file"],
        )

    @staticmethod
    def FUSION(a: str, b: str) -> str:
        """Either appends a with $ and b or just returns the one that isn't len 0"""
        if len(a) > 0 and len(b) > 0:
            sequence = a + "$" + b
        else:
            sequence = a + b
        return sequence

    @property
    def sequence_list(self) -> list[str]:
        return [f"{x},{self.abbreviation}" for x in self.instances.values()]

    @staticmethod
    def load_motif_json() -> list["MotifSequence"]:
        """load motifs.json file for curated motif collection"""
        json_path = Path(__file__).resolve().parent.joinpath("motifs.json")
        with open(json_path) as json_file:
            motif_list: list[MotifSequence] = json.load(
                json_file, object_hook=MotifSequence.from_dict
            )
        return motif_list

    def check_instance(self, annotation_string: str) -> bool:
        """Main regex searching function, checks the instance annotation string if it contains the instances associated motif. Removes any related/variation versions"""
        if len(annotation_string) < len(self.motif_name):
            return False
        searchd = re.search(
            pattern=self.motif_name.lower(),
            string=annotation_string.lower(),
        )
        if searchd is not None:
            return any(
                re.search(
                    pattern=x.lower(),
                    string=annotation_string.lower(),
                )
                for x in ["related", "mini"]
            )
        else:
            return False

    @property
    def instances(self) -> dict[str, str]:
        if hasattr(self, "_instances"):
            return self._instances
        else:
            return {}

    @instances.setter
    def instances(self, motif_instance: str) -> None:
        try:
            if hasattr(self, "_instances"):
                self._instances[motif_instance] = (
                    MotifSequence.get_rna3d_api_sequence(
                        motif_instance, self.loop_type
                    )
                )
            else:
                self._instances = {
                    motif_instance: MotifSequence.get_rna3d_api_sequence(
                        motif_instance, self.loop_type
                    )
                }
        except ValueError as e:
            logger.error(f"Error processing sequence {motif_instance}")
            logger.error(e)
        except ConnectionError as e:
            logger.critical("Could no reach RNA 3D Motif Atlas server")
            logger.error(e)
            raise

    @staticmethod
    def get_rna3d_api_sequence(instance_string: str, loop_type: str) -> str:
        """Main function for getting a sequence from the RNA 3D Motif Atlas API"""
        rna3datlas_api_call: str = (
            "http://rna.bgsu.edu/correspondence/pairwise_interactions_single?selection_type=loop_id&selection="
            + instance_string
        )
        api_return = call_api(rna3datlas_api_call).content.decode().split()
        nucleotides: list[str] = [x for x in api_return if "|" in x]
        Sequence = MotifSequence._extract_sequence_from_nucleotides(
            nucleotides, loop_type
        )
        if Sequence is not None:
            logger.debug(f"Retrieving {instance_string} was successful.")
            return Sequence
        if loop_type == "hairpin":
            raise ValueError(
                f"Error retreiving sequence for haripin motif instance {instance_string}, loop length was less than 3 nucleotides."
            )
        elif loop_type == "internal":
            raise ValueError(
                f"Error retreiving sequence for internal motif instance {instance_string}, sequence break was less than 5."
            )
        else:
            raise ValueError(
                f"Error retreiving sequence for haripin motif instance {instance_string}"
            )

    @staticmethod
    def _extract_sequence_from_nucleotides(
        nucleotide_list: list[str], loop_type: str
    ) -> str | None:
        """Subfunction for getting API sequences, extracts a sequence from a list of nucleotides"""
        match loop_type:
            case "hairpin":
                if len(nucleotide_list) - 2 >= 3:
                    return "".join(
                        MotifSequence.get_nucleotide_element(x, 3)
                        for x in nucleotide_list[1:-1]
                    )
                return None
            case "internal":
                sequence_break: int | None = MotifSequence.get_break(
                    nucleotide_list
                )
                if sequence_break is not None:
                    front: list[str] = nucleotide_list[1:sequence_break]
                    back: list[str] = nucleotide_list[
                        sequence_break + 2 : -1
                    ]
                    return MotifSequence.FUSION(
                        a="".join(
                            MotifSequence.get_nucleotide_element(x, 3)
                            for x in front
                        ),
                        b="".join(
                            MotifSequence.get_nucleotide_element(x, 3)
                            for x in back
                        ),
                    )
                else:
                    return None
            case _:
                logger.error(
                    "Unknown Looptype received for nucleotide extraction"
                )

    @staticmethod
    def get_break(nucleotide_list: list[str]) -> int | None:
        """Extracts sequence break from a list of nucleotide elements for Internal Loops. I cant use the break in the json file cause it does not account for bulges bases"""
        numbers: list[str] = [
            MotifSequence.get_nucleotide_element(x, 4)
            for x in nucleotide_list
        ]
        chains: list[str] = [
            MotifSequence.get_nucleotide_element(x, 2)
            for x in nucleotide_list
        ]
        for i in range(len(chains) - 1):
            if (
                chains[i + 1] != chains[i]
            ):  # Check if there is a sequence break with two different chain ids
                return i
            elif (
                abs(int(numbers[i + 1]) - int(numbers[i])) > 5
            ):  # check if there is a sequence break from
                return i
        return None

    @staticmethod
    def get_nucleotide_element(nucleotide: str, number: int) -> str:
        """Nucleotide element grabber that also converts non standard nucleotides using the nucleotide_conversion.json"""
        split: list[str] = nucleotide.split("|")
        element: str = split[number]
        # If you are grabbing the base from the nucleotide string
        if (number == 3) and element not in ["G", "C", "U", "A"]:
            element = MotifSequence.nucleotide_conversion(element)
        return element

    @staticmethod
    def nucleotide_conversion(nucleotide: str) -> str:
        """Conversion function from altered nucleotides back to their base nucleotides G/C/A/U"""
        if nucleotide in non_listed_conversions:
            return non_listed_conversions[nucleotide]
        with open(conversion_json_path, "r") as file:
            conversion_json = json.load(file)
        try:
            if conversion_json[nucleotide]["standard_base"][0] in [
                "G",
                "C",
                "A",
                "U",
            ]:
                return conversion_json[nucleotide]["standard_base"][0]
            else:
                raise ValueError(
                    f"Could not convert alternated nucleotide {nucleotide}, standard base is {conversion_json[nucleotide]['standard_base'][0]}"
                )
        except KeyError as e:
            raise ValueError(
                f"Could not convert alternated nucleotide {nucleotide}, raised {e}"
            )


@dataclass
class rna3d_motif:
    motif_id: str
    num_instances: int
    alignment: dict[str, list[str]]
    chainbreak: int
    common_name: str
    annotation: str
    bp_signature: str
    annotations: dict[str, str]

    @classmethod
    def from_dict(cls, rna3d_json_entry: dict[str, Any]) -> "rna3d_motif":
        """Main constructor for rna3d_motifs, used to create them from the hl/il.json files from the RNA 3D Motif Atlas"""
        return cls(
            motif_id=rna3d_json_entry["motif_id"],
            num_instances=rna3d_json_entry["num_instances"],
            alignment=rna3d_json_entry["alignment"],
            chainbreak=(
                rna3d_json_entry["chainbreak"][0]
                if rna3d_json_entry["chainbreak"]
                else 0
            ),
            common_name=rna3d_json_entry["common_name"],
            annotation=rna3d_json_entry["annotation"],
            bp_signature=rna3d_json_entry["bp_signature"],
            annotations=rna3d_json_entry["annotations"],
        )

    @staticmethod
    def load_rna3d_atlas(
        version: str = "current",
    ) -> list["rna3d_motif"]:
        """Retuns a list of all the current RNA 3D Motifs from the RNA 3D Motif Atlas as rna3d_motif objects"""
        version = version.replace("_", ".")
        current_hl_json = json.loads(
            call_api(
                f"http://rna.bgsu.edu/rna3dhub/motifs/release/hl/{version}/json"
            ).content.decode()
        )
        hairpin_objs: list[rna3d_motif] = [
            rna3d_motif.from_dict(hl_motif)
            for hl_motif in current_hl_json
            if "".join(hl_motif["annotations"].values())
        ]
        current_il_json = json.loads(
            call_api(
                f"http://rna.bgsu.edu/rna3dhub/motifs/release/il/{version}/json"
            ).content.decode()
        )
        internal_objs: list[rna3d_motif] = [
            rna3d_motif.from_dict(il_motif)
            for il_motif in current_il_json
            if "".join(il_motif["annotations"].values())
        ]
        return [*hairpin_objs, *internal_objs]


def assign_instances2MotifSequences(
    motif_sequences: list[MotifSequence],
    rna3d_mot_objs: list[rna3d_motif],
):
    """takes a list of MotifSequence objects and a list of rna3d_motif objects, iterates through the rna3d_motifs and assigns instances to corresponding MotifSequence objects"""
    # Split hairpins and internals to reduce the amount of necessary regex searches, not necessary but a decent enough time save for it to be worth a couple extra lines
    hairpins: list[MotifSequence] = [
        motif for motif in motif_sequences if motif.loop_type == "hairpin"
    ]
    internals: list[MotifSequence] = [
        motif for motif in motif_sequences if motif.loop_type == "internal"
    ]
    for motif in rna3d_mot_objs:
        if (
            motif.chainbreak == 0
        ):  # test case that is true for hairpin motifs
            regex_testing(motif, hairpins)
        else:  # test case that is true for internal motifs
            regex_testing(motif, internals)


def regex_testing(
    rna3d_mot: rna3d_motif, motif_sequences: list[MotifSequence]
):
    """Test function that takes an rna3d_motif and a list of MotifSequence objects and assigns annotated instances to the corresponding MotifSequence objects"""
    for instance in rna3d_mot.annotations:  # iterate instances
        if rna3d_mot.annotations[instance]:  # remove empty instances
            for motseq in motif_sequences:  # iterate motif_sequence objects
                if motseq.check_instance(rna3d_mot.annotations[instance]):
                    motseq.instances = instance
                    break


def call_api(call: str, attempts: int = 6) -> requests.Response:
    """Generic api call function, by default attempts 6 times to get a response from the API until it gives up"""
    for i in range(attempts):
        response = requests.get(call)
        if response.status_code == 200:
            logger.debug(f"API-request {call} successful")
            return response
        else:
            logger.error(
                f"Could not retrieve data in {i + 1} attempts. Request denied because: {response.reason}"
            )
    raise ConnectionError(
        f"Unable to retrieve RNA 3D Motif Atlas data from call: {call}"
    )


def sort_motif_sequences(
    list_of_MotifSequence_objs: list[MotifSequence],
) -> tuple[list[str], list[str], list[str]]:
    """Sorts MotifSequence objects into three separate lists for hairpins, internals and bulges (mainly done to separate internals and bulges)"""
    hairpins: list[str] = []
    internals: list[str] = []
    bulges: list[str] = []
    for motif in list_of_MotifSequence_objs:
        if motif.loop_type == "hairpin":
            hairpins.extend(motif.sequence_list)
        if motif.loop_type == "internal":
            for sequence in motif.sequence_list:
                if "$" in sequence:
                    internals.append(sequence)
                else:
                    bulges.append(sequence)
    return (hairpins, internals, bulges)


def write_csv(
    loop_type_sequences: list[str], loop_type: str, version: str
) -> None:
    """Writes all strings from a list of strings into a csv file named after the loop type"""
    setted = {x + "\n" for x in loop_type_sequences}
    mot_set_abbreviations = {x.split(",")[1].strip() for x in setted}
    output_path = (
        Path(__file__)
        .resolve()
        .parent.joinpath(
            "versions", "combined", f"{version}", f"rna3d_{loop_type}.csv"
        )
    )
    splitting_script = (
        Path(__file__).resolve().parent.joinpath("group_by_motif.sh")
    )
    Path.mkdir(
        self=Path(__file__)
        .resolve()
        .parent.joinpath("versions", "combined", f"{version}"),
        parents=False,
        exist_ok=True,
    )
    with open(
        file=output_path,
        mode="w+",
    ) as file:
        file.writelines(setted)
    motif_data = MotifSequence.load_motif_json()
    for mot in motif_data:
        if mot.abbreviation in mot_set_abbreviations:
            sep_output_path = (
                Path(__file__)
                .resolve()
                .parent.joinpath(
                    "versions",
                    "separated",
                    f"{version}_separated",
                    f"{loop_type}",
                    f"{mot.file}.csv",
                )
            )
            try:
                subprocess.run(
                    [
                        splitting_script,
                        output_path,
                        mot.abbreviation,
                        sep_output_path,
                    ],
                    check=True,
                )
            except ChildProcessError as e:
                logger.critical(
                    "Error during motif separation subprocess: " + str(e)
                )


def get_current_motif_version(attempts: int = 5) -> str:
    """retrieve current version of the RNA 3D Motif Atlas"""
    i = 0
    while i < attempts:
        response = requests.get(
            "http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json"
        )
        if response.status_code == 200:
            version = (
                response
                .headers["Content-disposition"]
                .split("=")[1]
                .split("_")[1][:-5]
                .strip()
            )
            logger.debug(
                f"Current RNA 3D Motif Atlas version is: {version}"
            )
            return version
        else:
            i += 1
    raise ConnectionError(
        f"Could not establish connection to API server in {attempts} attempts."
    )


def update_prep(version: str) -> bool:
    """Interactive update prep function parsing commandline for a specified version number and checking it against installed versions"""
    this_dir = Path(__file__).resolve().parent
    files_exist: list[bool] = [
        Path.is_file(
            this_dir.joinpath(
                "versions", "combined", version, "rna3d_bulges.csv"
            )
        ),
        Path.is_file(
            this_dir.joinpath(
                "versions", "combined", version, "rna3d_hairpins.csv"
            )
        ),
        Path.is_file(
            this_dir.joinpath(
                "versions", "combined", version, "rna3d_internals.csv"
            )
        ),
    ]
    if not all(files_exist):
        logger.debug("At least one of your motif files is missing")
        return True
    else:
        logger.debug(
            "All motif files of the requestsed version are available."
        )
        return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Python script for updating RNA 3D motifs",
        epilog="Who reads these anyways?",
    )
    parser.add_argument(
        "-v",
        "--version",
        type=str,
        help="Specify which RNA 3D Motif sequence version you want to download.",
        dest="version",
        default="current",
    )
    args = parser.parse_args()
    if args.version == "current":
        try:
            args.version = get_current_motif_version()
        except ConnectionError as error:
            logger.critical(error)
            raise
    args.version = args.version.replace(".", "_")
    return args


def update_hexdumps() -> None:
    filepath = Path(__file__).resolve().parent / "update_hexdump.sh"
    subprocess.run(
        [f"{filepath}"],
        cwd=Path(__file__).resolve().parent,
        check=True,
        shell=True,
    )


def compare_versions(
    requested_version: str, installed_version: str
) -> bool:
    return requested_version != installed_version


def currently_installed() -> str:
    try:
        with open(
            Path(__file__)
            .resolve()
            .parents[4]
            .joinpath("Extensions", "mot_header.hh"),
            "r",
        ) as file:
            header = file.readline()
    except FileNotFoundError:
        raise LookupError(
            "Unable to retrieve currently installed motif version"
        )
    version = re.search("[0-9].[0-9]+", header)
    if version is not None:
        return version.group()
    else:
        raise LookupError(
            "Unable to retrieve currently installed motif version"
        )


def uninteractive_update(requested_version: str) -> bool:  # type: ignore This function is for calling the updates from the main RNAmotiFold function so it is not used here
    """Updating function for RNA 3D Motif Sequence csv files, mainly for incorporation with other scripts (RNAmotiFold)"""
    logger.info(
        f"Uninteractive update process to version {requested_version} started"
    )
    if requested_version == "current":
        try:
            requested_version = get_current_motif_version()
        except ConnectionError as error:
            logger.critical(error)
            return False
    requested_version = requested_version.replace(".", "_")
    try:
        installed_version = currently_installed()
    except LookupError as e:
        logger.info(e)
        update_needed = True
    else:
        update_needed = compare_versions(
            requested_version, installed_version
        )
    if not update_needed:
        return False
    backd_up: bool = check_backups(requested_version=requested_version)
    if backd_up:
        update_hexdumps()
    else:
        updating: bool = update_prep(version=requested_version)
        if updating:
            return main(version_to_update_to=requested_version)
        else:
            return False
    return True


def check_backups(requested_version: str) -> bool:
    """Potential code for backup system if I ever have time to implement it (have a fully fledged system for using different RNA3D Motif Atlas Versions)."""
    versions_path: Path = (
        Path(__file__).resolve().parent.joinpath("versions", "combined")
    )
    for dir in get_dirs(versions_path):
        if dir == versions_path.joinpath(requested_version):
            logger.info(
                "Found requested version in backups, overwriting hexdump..."
            )
            return True
    logger.info(
        "Could not find a backup with that version, attempting to fetch motifs from RNA 3D Motif Atlas..."
    )
    return False


def get_dirs(root: str | Path) -> Generator[Path, None, None]:
    """
    Generator, which recursively yields all subdirectories of a directory
    """
    for path in Path(root).rglob("*"):
        if path.is_dir():
            yield path


def get_files(root: str | Path) -> Generator[Path, None, None]:
    """Generator the recursively yields all csv files of a directory"""
    for path in Path(root).rglob("*"):
        if path.is_file() and path.suffix == ".csv":
            yield path


def main(version_to_update_to: str) -> bool:
    """Main body of the script, gets the version to update to either interactively from the user via the command line or uninteractively just updates to the current RNA 3D Motif Atlas version"""
    try:
        motif_list: list[MotifSequence] = MotifSequence.load_motif_json()
        logger.info("Retrieved motif json")
        rna3d_motif_jsons: list[rna3d_motif] = rna3d_motif.load_rna3d_atlas(
            version_to_update_to
        )
        logger.info(
            f"Retrieved RNA 3D Motif Atlas version {version_to_update_to}"
        )
        assign_instances2MotifSequences(
            motif_sequences=motif_list, rna3d_mot_objs=rna3d_motif_jsons
        )
        hairpins, internals, bulges = sort_motif_sequences(
            list_of_MotifSequence_objs=motif_list
        )
        write_csv(
            loop_type_sequences=hairpins,
            loop_type="hairpins",
            version=version_to_update_to,
        )
        write_csv(
            loop_type_sequences=internals,
            loop_type="internals",
            version=version_to_update_to,
        )
        write_csv(
            loop_type_sequences=bulges,
            loop_type="bulges",
            version=version_to_update_to,
        )
        update_hexdumps()
        return True
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        return False


if __name__ == "__main__":
    args = parse_args()
    main(version_to_update_to=args.version)
