import os
import subprocess
import tempfile

from pichemist.config import ACD_METHOD
from pichemist.config import PKA_LIMITS


class ACDPKaException(Exception):
    pass


class ACDPKaCalculator(object):
    """Uses ACD perceptabat to calculate pKa values."""

    def __init__(self):
        self.input_filepath = self._get_temp_filepath(".smi")
        self.output_filepath = self._get_temp_filepath(".out")
        self.pka_flag = self._get_pka_flag()

    def _get_pka_flag(self):
        return ACD_METHOD.value

    def _get_temp_filepath(self, suffix):
        """Gets a temporary filepath and closes the file."""
        file = tempfile.NamedTemporaryFile(suffix=suffix)
        path = file.name
        file.close()
        return path

    def _prepare_temp_input_file(self, smi_list):
        """Prepares the input file with SMILES."""
        with open(self.input_filepath, "w") as f:
            i = 0
            for smi in smi_list:
                i += 1
                f.write(f"{smi} tmpname{i}\n")

    def _get_temp_output_filepath(self):
        """Gets a temporary filepath and closes the file."""
        file = tempfile.NamedTemporaryFile(suffix=".out")
        path = file.name
        file.close()
        return path

    def _build_command(self):
        """Builds the perceptabat command to run pKa prediction."""
        return [
            "perceptabat",
            f"-TFNAME{self.output_filepath}",
            self.pka_flag,
            "-TPKA",
            self.input_filepath,
        ]

    def _get_status_output(self, *args, **kwargs):
        p = subprocess.Popen(*args, **kwargs)
        stdout, stderr = p.communicate()
        return p.returncode, stdout, stderr

    def _run_acd_exe(self, cmd):
        status, stdout, stderr = self._get_status_output(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        if status != 0:
            raise ACDPKaException(
                f"""Error while running the ACD subprocess:
                                  status: {status}
                                  command: {cmd}
                                  stdout: {stdout}
                                  stderr: {stderr}"""
            )

    def _check_existence_output(self):
        if not os.path.isfile(self.output_filepath):
            raise ACDPKaException("ACD did not generate the output file.")

    def _delete_temp_files(self):
        """Deletes the temp files of the instance."""
        for path in [self.input_filepath, self.output_filepath]:
            if os.path.exists(path):
                os.remove(path)

    def _parse_output(self, smi_list):
        """Parses the output of the software."""
        with open(self.output_filepath, "r") as f:
            base_pkas = list()
            acid_pkas = list()
            f.readline()  # skip first line
            for line in f.readlines():
                ln = line.split()
                mol_idx = int(ln[0])
                if "ACD_pKa_Apparent" in ln[1]:
                    pka = float(ln[2])
                if "ACD_pKa_DissType_Apparent" in ln[1]:
                    if ln[2] in ["MB", "B"]:
                        if pka > PKA_LIMITS["base_1"] and pka < PKA_LIMITS["base_2"]:
                            base_pkas.append((pka, smi_list[mol_idx - 1]))
                    if ln[2] in ["MA", "A"]:
                        if pka > PKA_LIMITS["acid_1"] and pka < PKA_LIMITS["acid_2"]:
                            acid_pkas.append((pka, smi_list[mol_idx - 1]))
        return (base_pkas, acid_pkas)

    def calculate_pka_from_list(self, smi_list):
        """Calculates the pKa values of a list of SMILES."""
        # Run the calculations
        self._prepare_temp_input_file(smi_list)
        cmd = self._build_command()
        self._run_acd_exe(cmd)
        self._check_existence_output()
        results = self._parse_output(smi_list)
        self._delete_temp_files()
        return results
