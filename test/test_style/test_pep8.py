import os
import unittest

import pycodestyle


class TestPep8(unittest.TestCase):
    """Run PEP8 on all files in this directory and subdirectories."""

    _exclude_dirs = ['subrepos', 'utility_scripts', '.venv', '.pytest_cache', 'pheconstructors']

    def setUp(self):
        self.base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.exclude_dirs = [os.path.join(self.base_dir, d) for d in self._exclude_dirs]

    def test_pep8(self):
        # This is here so that we can search for the start of errors in the logs
        print("Checking for PEP8 errors")
        style = pycodestyle.StyleGuide(max_line_length=120)
        errors = 0
        for root, directory, files in os.walk(self.base_dir):
            directory[:] = [d for d in directory if os.path.join(root, d) not in self.exclude_dirs]
            python_files = [f for f in files if f.endswith('.py')]
            for pf in python_files:
                check = style.check_files([os.path.join(root, pf)])
                errors += check.file_errors
        self.assertEqual(errors, 0)


if __name__ == "__main__":
    unittest.main()
