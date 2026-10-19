import unittest

import RNAmotiFold.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_rna3d_motifs as motifs
from requests import Response

##How To Testing:
# ARRANGE -> ACT -> ASSERT
# Arrange a test case, class in a specific state,
# Act out a function that the object might do in the sate
# Assert that the right thing was returned

# Test behavior, not implementation details:
# Set tests up in a way that even if implementation of a function changes
# the test should still stand as long as the result doesn't change
#
# Each test should be for one specific behavior!
#
#
# Test edge cases and failures!


class TestBackUpChecking(unittest.TestCase):
    def test_backup_check(self):
        """Test for backup checking"""
        # Arrange
        requested_version = "4_2"
        # Act
        returnval = motifs.check_backups(requested_version)
        # Assert
        self.assertIs(returnval, True)

    def test_backup_check_fail(self):
        """Test for non back'd up version"""
        # Arrange
        request_version = "3_8"
        # Act
        returnval = motifs.check_backups(requested_version=request_version)
        # Assert
        self.assertIs(returnval, False)


class TestServerAvailable(unittest.TestCase):
    def test_RNA3DMotifAtlasAvailable(self):
        """Tests if RNA 3D Motif Atlas Server is reachable"""
        # Arrange
        api_call = (
            f"http://rna.bgsu.edu/rna3dhub/motifs/release/il/4.11/json"
        )
        # Act
        returnval: Response = motifs.call_api(api_call, 5)
        # Assert
        self.assertIsInstance(returnval, Response)


class TestNucleotideProcessing(unittest.TestCase):
    def test_getting_rna3d_api_sequence(self):
        """Test if we can get and process a motif from the server"""
        # Arrange
        instance = "IL_3V7E_006"
        type = "internal"
        # Act
        sequence = motifs.MotifSequence.get_rna3d_api_sequence(
            instance, type
        )
        # Assert
        self.assertEqual(sequence, "GGA$GAUGAA")

    def test_sequence_extraction_internal(self):
        """Test only the sequence extraction from the motif processing pipeline for internals"""
        # Arrange
        # Nucleotide Instance: PDB ID | MODEL NUMBER | CHAIN ID | NUCLEOTIDE | POS | Atom Name | Alternate ID | Insertion Code | Symmetry Operation
        nucleotides = [
            "Test|1|A|G|2",
            "Test|1|A|C|3",
            "Test|1|A|C|4|||C",
            "Test|1|A|G|5",
            "Test|1|A|C|11",
            "Test|1|A|G|12",
            "Test|1|A|G|11",
            "Test|1|A|C|13|||G",
        ]
        loop_type = "internal"
        # Act
        sequence = motifs.MotifSequence._extract_sequence_from_nucleotides(
            nucleotides, loop_type
        )
        self.assertEqual(sequence, "CC$GG")

    def test_sequence_extraction_hairpin(self):
        """Test only the sequence extraction from the motif processing pipeline for hairpins"""
        # Arrange
        nucleotides = [
            "Test|1|A|G|1",
            "Test|1|A|C|2",
            "Test|1|A|G|3|||C",
            "Test|1|A|A|4",
            "Test|1|A|G|5",
            "Test|1|A|A|6",
            "Test|1|A|G|4",
            "Test|1|A|C|4",
        ]
        loop_type = "hairpin"
        # Act
        sequence = motifs.MotifSequence._extract_sequence_from_nucleotides(
            nucleotides, loop_type
        )
        # Assert
        self.assertEqual(sequence, "CGAGAG")

    def test_break_finder(self):
        """Test break finder for internal loops with legit internal loop"""
        # Arrange
        nucleotides = [
            "Test|1|A|G|1",
            "Test|1|A|C|2",
            "Test|1|A|C|3|||C",
            "Test|1|A|G|4",
            "Test|1|A|C|10",
            "Test|1|A|G|11",
            "Test|1|A|G|12",
            "Test|1|A|C|12|||G",
        ]
        # Act
        seq_break = motifs.MotifSequence.get_break(nucleotides)
        # Assert
        self.assertEqual(seq_break, 3)

    def test_break_finder_nobreak(self):
        """Test break finder for internal loops with faulty sequence"""
        # Arrange
        nucleotides = [
            "Test|1|A|G|1",
            "Test|1|A|C|2",
            "Test|1|A|C|3|||C",
            "Test|1|A|G|4",
            "Test|1|A|C|7",
            "Test|1|A|G|10",
            "Test|1|A|G|11",
            "Test|1|A|C|12|||G",
        ]
        # Act
        seq_break = motifs.MotifSequence.get_break(nucleotides)
        # Assert
        self.assertIsNone(seq_break)


if __name__ == "__main__":
    unittest.main()
