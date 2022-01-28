import argparse
import unittest

from ukb_pret.cli import preprocess_command_line_inputs


class TestCommandLineInterface(unittest.TestCase):

    def test_no_prs_files(self):
        input = ['--pheno-file', 'a/b/c.csv']
        with self.assertRaises(SystemExit) as se:
            preprocess_command_line_inputs(input)

    def test_no_pheno_file(self):
        input = ['--prs-files', 'a/b/c.csv']
        with self.assertRaises(SystemExit) as se:
            preprocess_command_line_inputs(input)

    def test_multiple_prs_files(self):
        input = ['--prs-files', 'a/b/c.csv', '/d/e/f.csv', '--pheno-file', 'test.pheno']
        preprocess_command_line_inputs(input)
