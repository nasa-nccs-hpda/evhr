import sys

import unittest
from unittest.mock import MagicMock
from unittest.mock import Mock

sys.modules['ipyfilechooser'] = Mock()
sys.modules['ipyleaflet'] = Mock()
sys.modules['ipywidgets'] = Mock()
sys.modules['IPython.display'] = Mock()
sys.modules['localtileserver'] = Mock()

from evhr.view.notebook.EvhrToaDashboard import EvhrToaDashboard


class EvhrDashboardTestCase(unittest.TestCase):

    # -------------------------------------------------------------------------
    # test_check_number_of_files_matching_regex
    # -------------------------------------------------------------------------
    def test_check_number_of_files_matching_regex(self):
        # Create an instance of your class
        dashboard = EvhrToaDashboard()

        # Define the file list and catalog ID for testing
        file_list = ['/path/to/file1/12345',
                     '/path/to/file2/12345',
                     '/path/to/file3/12345']

        cat_id = '12345'

        # Call the method to be tested
        dashboard._check_number_of_files_matching_regex(file_list, cat_id)

        # Assert that the warning value is set correctly
        expected_warning = f'Found {len(file_list)} TOAs' + \
            f' matching catelog ID {cat_id}'
        self.assertEqual(dashboard._warning.value, expected_warning)

        # Define a different file list with fewer files
        file_list = []

        # Call the method again with the new file list
        dashboard._check_number_of_files_matching_regex(file_list, cat_id)

        # Assert that the error message is set correctly
        expected_error = '<b>ERROR: Could not find any TOAs' + \
            f' matching catelog ID {cat_id}</b>'

        self.assertEqual(dashboard._warning.value, expected_error)

    # -------------------------------------------------------------------------
    # test_get_tag_not_found
    # -------------------------------------------------------------------------
    def test_get_tag_not_found(self):
        # Create an instance of your class
        dashboard = EvhrToaDashboard()

        # Define the tag string and XML file name for testing
        tag_str = 'non_existing_tag'
        xml_file_name = 'example.xml'

        # Create a mock XML tree
        xml_tree = MagicMock()

        # Mock the getroot method to return a mock element
        root_element = MagicMock()
        xml_tree.getroot.return_value = root_element

        # Mock the find method to return None
        root_element.find.return_value = None

        # Mock the _warning.value attribute
        dashboard._warning.value = ''

        # Call the method to be tested and assert that it raises a RuntimeError
        with self.assertRaises(RuntimeError):
            dashboard._get_tag(tag_str, xml_tree, xml_file_name)

        # ---
        # Assert that the _warning.value attribute is set
        # with the expected error message
        # ---
        expected_error_msg = \
            f'<b>Unable to locate the "{tag_str}" tag in {xml_file_name}</b>'
        self.assertEqual(dashboard._warning.value, expected_error_msg)

        # Assert that the getroot method was called once
        xml_tree.getroot.assert_called_once()

        # ---
        # Assert that the find method was called once with the
        # expected arguments
        # ---
        root_element.find.assert_called_once_with(tag_str)

    # -------------------------------------------------------------------------
    # test_calculate_expected_time
    # -------------------------------------------------------------------------
    def test_calculate_expected_time(self):
        # Create an instance of your class
        dashboard = EvhrToaDashboard()

        # Define the image size for testing
        img_size = 10_000 * 10_000

        # Define the expected time coefficient
        expected_time_coeff = 0.00000707

        # Set the EVHR_EXPECTED_TIME_COEFF to the expected time coefficient
        dashboard.EVHR_EXPECTED_TIME_COEFF = expected_time_coeff

        # Call the method to be tested
        result = dashboard._calculate_expected_time(img_size)

        # ---
        # Calculate the expected time based on the image size and
        # expected time coefficient
        # ---
        expected_time = img_size * expected_time_coeff

        # Assert that the result is the expected time
        self.assertEqual(result, expected_time)
