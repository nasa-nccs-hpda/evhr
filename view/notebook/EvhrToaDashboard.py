import concurrent.futures
import functools
import glob
import os
import pathlib
import subprocess
import time
from typing import Tuple
import xml.etree.ElementTree as ET

from ipyfilechooser import FileChooser
from ipyleaflet import Map, basemaps, ScaleControl, LayersControl
from ipyleaflet import FullScreenControl
import ipywidgets as widgets
from IPython.display import display, update_display
from localtileserver import TileClient, get_leaflet_tile_layer

from core.model.NotebookImageHelper import NotebookImageHelper


SINGULARITY_DIR: str = '/explore/nobackup/people/iluser/' + \
    'ilab_containers/prod_sandboxes'

EVHR_CONTAINER_NAME_REGEX: str = '*evhr*'

SING_EXEC_PRE_STR: str = 'singularity exec -B /explore,/panfs,/css,/nfs4m,'

EVHR_TOA_CLV_EXEC: str = 'python /usr/local/ilab/evhr/view/evhrToaCLV.py'

SOURCE_BASH_RC_CMD: str = 'source ~/.bashrc'

MODULE_LOAD_CMD: str = 'module load singularity'

EVHR_TOA_CLV_ARGS: dict = {"output_dir": '-o',
                           "pan_res": '--pan_res',
                           "pan_sharpen": '--pan_sharpen',
                           "scene": '--scenes'}

CATELOG_ID_TAG: str = './IMAGE/CATID'

IMD: str = 'IMD'

TIL: str = 'TIL'

WIDTH_TIL_TAG: str = 'TILESIZEX'

HEIGHT_TIL_TAG: str = 'TILESIZEY'

ACCEPTABLE_FILE_COUNT: int = 1

RED_BAND_ID: int = 3

GREEN_BAND_ID: int = 2

BLUE_BAND_ID: int = 1

EVHR_NODATA_VAL: int = -10_001

EVHR_EXPECTED_TIME_COEFF: float = 0.00000895


# -------------------------------------------------------------------------
# EvhrToaDashboard
# -------------------------------------------------------------------------
class EvhrToaDashboard(object):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self) -> None:

        self._title = widgets.HTML(
            value="<h1><b>EVHR TOA Dashboard</b><h1>")

        self._input_file_chooser = FileChooser(
            title='<b>EVHR Input</b>',
            layout=widgets.Layout(width='auto'))

        self._pan_res_drop_down = widgets.Dropdown(
            options=['0.5', '1'],
            value='0.5',
            description='Pan Resolution:',
            disabled=False)
        self._pan_sharpen_checkbox = widgets.Checkbox(
            description='Pan Sharpen')

        self._output_dir_chooser = FileChooser(title='<b>EVHR Output Dir</b>',
                                               layout=widgets.Layout(
                                                   width='auto'),
                                               show_only_dirs=False)

        self._submit_button = widgets.Button(description='Submit')

        self._warning = widgets.HTML(value="")

        self._imd_tag = None

        self._til_tag = None

        self._cla_widgets = None

        self._display_id = str(time.time())

    # -------------------------------------------------------------------------
    # render
    # -------------------------------------------------------------------------
    def render(self):
        """Renders the EVHR dashboard. A dashboard which is used to run the
        EVHR TOA program through a GUI in a Jupyter notebook.
        """

        self._cla_widgets = self._organize_widgets()

        self._submit_button.on_click(self._launch_submit_and_display_process)

        display(self._cla_widgets, display_id=self._display_id)

    # -------------------------------------------------------------------------
    # _organize_widgets
    # -------------------------------------------------------------------------
    def _organize_widgets(self) -> widgets.widget_box.VBox:
        """Organizes the widgets that make up the dashboard

        Returns:
            widgets.widget_box.VBox: widgets box which contains all widgets
            instantiated
        """

        running_args_widgets_list = [self._pan_res_drop_down,
                                     self._pan_sharpen_checkbox]

        running_args_box = widgets.HBox(running_args_widgets_list)

        warning_widgets_list = [self._warning]

        warning_box = widgets.HBox(warning_widgets_list, box_style='warning')

        cla_widgets_list = [self._title,
                            self._input_file_chooser,
                            running_args_box,
                            self._output_dir_chooser,
                            warning_box,
                            self._submit_button]

        cla_widgets_box = widgets.VBox(cla_widgets_list)

        return cla_widgets_box

    # -------------------------------------------------------------------------
    # _launch_submit_and_display_process
    # -------------------------------------------------------------------------
    def _launch_submit_and_display_process(self, button: widgets.Button):
        """When the submit button is clicked, this function is run which starts
        the run process

        Args:
            button (widgets.Button): Button which is tied to the on-click
            function
        """

        input_file_path = pathlib.Path(self._input_file_chooser.value)

        if not input_file_path.exists():

            post_msg_str = 'does not exist. Choose a proper file.</b>'

            self._warning.value = f'<B> {input_file_path} {post_msg_str}'

        cla_dict = {
            'scene': str(input_file_path),
            'pan_res': self._pan_res_drop_down.value,
            'pan_sharpen': self._pan_sharpen_checkbox.value,
            'output_dir': pathlib.Path(self._output_dir_chooser.value),
        }

        self._submit_button.description = 'Submitting'

        evhr_toa_file = self._evhr_launch_process(cla_dict)

        self._submit_button.description = 'Done'

        self._render_update_map_display(evhr_toa_image_path=evhr_toa_file)

    # -------------------------------------------------------------------------
    # _evhr_launch_process
    # -------------------------------------------------------------------------
    def _evhr_launch_process(self, cla_dict: dict) -> str:
        """Given a command-line dict, sets up and launches EVHR command

        Args:
            cla_dict (dict): Command-line argument dictionary

        Returns:
            str: EVHR TOA image file path
        """

        self._get_til_and_imd_tag_xml(cla_dict['scene'])

        evhr_toa_command_str = self._make_singularity_command(cla_dict)

        pbar, expected_time = self._setup_progress_bar()

        self._render_update_progress_bar(pbar)

        @EvhrToaDashboard.progress_bar(float_progress=pbar,
                                       expected_time=expected_time)
        def launch_subprocess(command: str):
            self._launch_subprocess(command)

        launch_subprocess(evhr_toa_command_str)

        evhr_toa_file = self._find_output_file(cla_dict['output_dir'])

        return evhr_toa_file

    # -------------------------------------------------------------------------
    # _get_xml_tree
    # -------------------------------------------------------------------------
    def _get_til_and_imd_tag_xml(self, file_name: str) -> None:
        """Opens and reads xml file, searches for IMD and TIL tags in the xml
        tree.

        Args:
            file_name (str): XML file name

        Raises:
            RuntimeError: If companion file is not TIF or NTF, raises error.
            FileNotFoundError: If XML file does not exist, raises error.
        """

        extension = os.path.splitext(file_name)[1]

        if extension != '.ntf' and extension != '.tif':
            raise RuntimeError('{} is not a NITF or TIFF'.format(file_name))

        # Ensure the XML file exists.
        xml_file_name = file_name.replace(extension, '.xml')

        if not os.path.exists(xml_file_name):
            error_msg = '{} does not exist'.format(xml_file_name)
            self._warning.value = f"<b>{error_msg}</b>"
            raise FileNotFoundError()

        # Some data members require the XML file counterpart to the TIF.
        xml_tree = ET.parse(xml_file_name)

        # Get imd tag for getting the cat ID
        self._imd_tag = self._get_tag(IMD, xml_tree, xml_file_name)

        # Get til tag for generating progress bar
        self._til_tag = self._get_tag(TIL, xml_tree, xml_file_name)

    # -------------------------------------------------------------------------
    # _get_tag
    # -------------------------------------------------------------------------
    def _get_tag(self,
                 tag_str: str,
                 xml_tree: ET,
                 xml_file_name: str):
        """Searches for a tag in an xml tree

        Args:
            tag_str (str): tag string to search for
            xml_tree (ET): xml tree to search
            xml_file_name (str): file name of the xml file

        Raises:
            RuntimeError: Raises if the tag is not found in the xml tree.

        Returns:
            _type_: xml sub-tree matching tag
        """

        xml_tag = xml_tree.getroot().find(tag_str)

        if xml_tag is None:
            error_msg = \
                f'Unable to locate the "{tag_str}" tag in ' + xml_file_name
            self._warning.value = f"<b>{error_msg}</b>"
            raise RuntimeError()

        return xml_tag

    # -------------------------------------------------------------------------
    # _make_singularity_command
    # -------------------------------------------------------------------------
    def _make_singularity_command(self, cla_dict: dict) -> str:
        """Builds the entire EVHR invocation demand including sourcing
        `bashrc` and singularity execs.

        Args:
            cla_dict (dict): dict containing EVHR TOA command-line arguments

        Returns:
            str: EVHR TOA invocation command
        """

        evhr_container_path = self._get_evhr_container_path()

        cmd = f"{SOURCE_BASH_RC_CMD}"

        cmd = f"{cmd} && {MODULE_LOAD_CMD} && {SING_EXEC_PRE_STR}"

        cmd = f"{cmd} {evhr_container_path} {EVHR_TOA_CLV_EXEC}"

        for cla_key in cla_dict:

            argument_value = cla_dict[cla_key]

            argument_key = EVHR_TOA_CLV_ARGS[cla_key]

            if not argument_value:
                continue

            if type(argument_value) == bool:
                cmd = f"{cmd} {argument_key}"

            else:
                cmd = f"{cmd} {argument_key} {argument_value}"

        cmd = f"({cmd})"

        return cmd

    # -------------------------------------------------------------------------
    # _get_evhr_container_path
    # -------------------------------------------------------------------------
    def _get_evhr_container_path(self) -> pathlib.Path:
        """Find the container which matches the container regex. Returns the
        path to the container.

        Returns:
            pathlib.Path: evhr container path
        """

        evhr_container_regex = os.path.join(
            SINGULARITY_DIR, EVHR_CONTAINER_NAME_REGEX)

        glob_attempt = glob.glob(evhr_container_regex)

        self._check_glob(glob_attempt, evhr_container_regex)

        evhr_container = pathlib.Path(glob_attempt[0])

        return evhr_container

    # -------------------------------------------------------------------------
    # _check_glob
    # -------------------------------------------------------------------------
    def _check_glob(self, glob_attempt: list, evhr_regex: str) -> None:
        """Check if any container paths match the regex. Raises runtime error
        if not.

        Args:
            glob_attempt (list): List of files matching description
            evhr_regex (str): Regex to find files matching to

        Raises:
            FileNotFoundError: Raised if no files matched the regex.
        """

        if len(glob_attempt) < 1:

            error_msg = f'Could not find container at {evhr_regex}'
            self._warning.value = f"<b>{error_msg}</b>"
            raise FileNotFoundError(error_msg)

    # -------------------------------------------------------------------------
    # getCatalogId()
    # -------------------------------------------------------------------------
    def _get_catalog_id(self) -> str:

        return self._imd_tag.findall(CATELOG_ID_TAG)[0].text

    # -------------------------------------------------------------------------
    # _setup_progress_bar()
    # -------------------------------------------------------------------------
    def _setup_progress_bar(self) -> Tuple:
        """Get image size and expected EVHR runtime, instantiate a widgets
        progress bar. Returns progres bar and expected time.

        Returns:
            Tuple: (progress_bar, expected_time)
        """

        image_size = self._get_image_size()

        expected_time = self._calculate_expected_time(image_size)

        expected_time = int(expected_time)

        progress = widgets.FloatProgress(value=0.0, min=0.0, max=expected_time)

        return progress, expected_time

    # -------------------------------------------------------------------------
    # _get_image_size()
    # -------------------------------------------------------------------------
    def _get_image_size(self) -> int:
        """Given the XML tag with image width and height, calculate total area
        in pixel. Returns the image total area (image size)

        Returns:
            int: image size
        """

        # Search for the specified tags in the XML file
        image_width = self._til_tag.findall(WIDTH_TIL_TAG)[0].text

        image_height = self._til_tag.findall(HEIGHT_TIL_TAG)[0].text

        image_size = int(image_width) * int(image_height)

        return image_size

    # -------------------------------------------------------------------------
    # _calculate_expected_time()
    # -------------------------------------------------------------------------
    def _calculate_expected_time(self, img_size: int) -> float:
        """Given an image size, calculate an estimate (in seconds) how long
        EVHR is expected to run.

        Args:
            img_size (int): TOA image size

        Returns:
            float: expected time (in seconds)
        """

        expected_time = img_size * EVHR_EXPECTED_TIME_COEFF

        return expected_time

    # -------------------------------------------------------------------------
    # _render_update_progress_bar()
    # -------------------------------------------------------------------------
    def _render_update_progress_bar(self, pbar):
        """Update the application to display the progress bar.

        Args:
            pbar (_type_): progress bar
        """

        updated_display_vbox = widgets.VBox([self._cla_widgets, pbar])

        update_display(updated_display_vbox, display_id=self._display_id)

    # -------------------------------------------------------------------------
    # _launch_subprocess
    # -------------------------------------------------------------------------
    @staticmethod
    def _launch_subprocess(command: str) -> Tuple[str, str]:
        """Launch a subprocess through subprocess module.

        Args:
            command (str): command string to launch

        Raises:
            RuntimeError: Executing of subprocess failed
            RuntimeError: An error occuring while launching a subprocess

        Returns:
            Tuple[str, str]: output, error strings
        """

        try:

            process = subprocess.Popen(command,
                                       shell=True,
                                       stderr=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       close_fds=True)

            output, error = process.communicate()

            return output, error

        except subprocess.CalledProcessError as e:

            error_msg = f"Subprocess failed. Command: {command}. Error: {e}"

            raise RuntimeError(error_msg)

        except Exception as e:

            error_msg = "An error occurred while launching the " + \
                f"subprocess. Command: {command}. Error: {e}"

            raise RuntimeError(error_msg)

    # -------------------------------------------------------------------------
    # _find_output_file()
    # -------------------------------------------------------------------------
    def _find_output_file(self, output_dir: pathlib.Path) -> pathlib.Path:
        """Given a catelog ID string, find the corresponding TOA file path
        that contains the catalog ID.

        Args:
            output_dir (pathlib.Path): output directory

        Returns:
            pathlib.Path: TOA file
        """

        cat_id = self._get_catalog_id()

        regex = f'5-toas/*{cat_id}*-toa.tif'

        files_matching_regex = list(output_dir.glob(regex))

        self._check_number_of_files_matching_regex(
            files_matching_regex, cat_id)

        file_matching_regex_to_return = pathlib.Path(files_matching_regex[0])

        return file_matching_regex_to_return

    # -------------------------------------------------------------------------
    # _check_number_of_files_matching_regex
    # -------------------------------------------------------------------------
    def _check_number_of_files_matching_regex(self, file_list: list,
                                              cat_id: str) -> None:
        """Checks if more than one TOA file exists that contains the catalog
        ID string. Warns if more than one file exists. Warns if no files match.

        Args:
            file_list (list): list of paths
            cat_id (str): catalog ID
        """

        number_of_files = len(file_list)

        if number_of_files > ACCEPTABLE_FILE_COUNT:

            warning_msg = f'Found {number_of_files} TOAs' + \
                f' matching catelog ID {cat_id}'

            self._warning.value = warning_msg

        if number_of_files < ACCEPTABLE_FILE_COUNT:

            error_msg = '<b>ERROR: Could not find any TOAs' + \
                f' matching catelog ID {cat_id}</b>'

            self._warning.value = error_msg

    # -------------------------------------------------------------------------
    # _render_update_map_display
    # -------------------------------------------------------------------------
    def _render_update_map_display(self,
                                   evhr_toa_image_path: pathlib.Path) -> None:
        """Renders a ipyleaflet map with the TOA image.

        Args:
            evhr_toa_image_path (pathlib.Path): Path to toa image.
        """

        evhr_toa_image_helper = self._init_image_helper(evhr_toa_image_path)

        evhr_toa_image_layer, evhr_toa_image_client = self._init_tile(
            evhr_toa_image_path,
            evhr_toa_image_helper)

        m = Map(
            center=evhr_toa_image_client.center(),
            zoom=evhr_toa_image_client.default_zoom,
            basemap=basemaps.Esri.WorldImagery,
            scroll_wheel_zoom=True,
            keyboard=True,
            dragging=True,
            layout=widgets.Layout(height='600px')
        )

        m.add_layer(evhr_toa_image_layer)
        m.add_control(ScaleControl(position='bottomleft'))
        m.add_control(LayersControl(position='topright'))
        m.add_control(FullScreenControl())

        updated_display_vbox = widgets.VBox([self._cla_widgets, m])

        update_display(updated_display_vbox, display_id=self._display_id)

    # -------------------------------------------------------------------------
    # _init_image_helper
    # -------------------------------------------------------------------------
    def _init_image_helper(
            self,
            evhr_toa_image_path: pathlib.Path) -> NotebookImageHelper:
        """Initialize image helper with TOA properties

        Args:
            evhr_toa_image_path (pathlib.Path): Path to TOA image.

        Returns:
            NotebookImageHelper: NotebookImageHelper for TOA
        """

        evhr_toa_image_helper = NotebookImageHelper()

        evhr_toa_image_helper.initFromFile(inputFile=evhr_toa_image_path,
                                           noDataValue=EVHR_NODATA_VAL,
                                           redBandId=RED_BAND_ID,
                                           greenBandId=GREEN_BAND_ID,
                                           blueBandId=BLUE_BAND_ID)

        return evhr_toa_image_helper

    # -------------------------------------------------------------------------
    # _init_tile
    # -------------------------------------------------------------------------
    def _init_tile(self,
                   evhr_toa_image_path: pathlib.Path,
                   evhr_toa_image_helper: NotebookImageHelper) -> Tuple:
        """Initialize tile clients and layers for TOA image.

        Args:
            evhr_toa_image_path (pathlib.Path): _description_
            evhr_toa_image_helper (NotebookImageHelper): _description_

        Returns:
            Tuple: leaflet_tile_layer, leaflet_tile_client
        """

        evhr_toa_image_tile_client = TileClient(evhr_toa_image_path)

        evhr_toa_layer = get_leaflet_tile_layer(
            evhr_toa_image_tile_client,
            nodata=evhr_toa_image_helper._noDataValue,
            vmin=evhr_toa_image_helper._minValue,
            vmax=evhr_toa_image_helper._maxValue,
            name=evhr_toa_image_path.stem,
            band=evhr_toa_image_helper.getRgbIndexList()
        )

        return evhr_toa_layer, evhr_toa_image_tile_client

    # -------------------------------------------------------------------------
    # progress_bar
    # -------------------------------------------------------------------------
    @staticmethod
    def progress_bar(float_progress: widgets.widget_float.FloatProgress,
                     expected_time: int,
                     increments: int = 100) -> None:
        """Progress bar function that given a ipywidgets progress bar, will
        update the progress bar asynchronously while the function is executed.

        Uses expected time and increments to determine how often to update the
        progress bar.

        From: https://stackoverflow.com/questions/59013308/\
            python-progress-bar-for-non-loop-function

        Args:
            float_progress (widgets.widget_float.FloatProgress): widgets float
                progress bar.
            expected_time (int): amount of time (in seconds) the wrapped
                function is expected to take
            increments (int, optional): In what increment to update the
                progress bar. Defaults to 100.

        Example:
            ```python
            progress = widgets.FloatProgress(value=0.0, min=0.0, max=10)

            @progress_bar(float_progress=progress, expected_time=10)
            def test_func():
                time.sleep(8)
                return "result"

            display(progress)
            print(test_func())  # prints "result"
            ```
        """

        def _progress_bar(func):

            def timed_progress_bar(future,
                                   float_progress,
                                   expected_time,
                                   increments=100):
                """
                Display progress bar for expected_time seconds.
                Complete early if future completes.
                Wait for future if it doesn't complete in expected_time.
                """
                interval = expected_time / increments
                for i in range(increments - 1):
                    if future.done():
                        # finish the progress bar
                        # not sure if there's a cleaner way to do this?
                        float_progress.value = expected_time
                        float_progress.bar_style = 'success'
                        return
                    else:
                        time.sleep(interval)
                        float_progress.value = float_progress.value + interval
                # if the future still hasn't completed, wait for it.
                future.result()
                float_progress.value = expected_time
                float_progress.bar_style = 'success'

            @functools.wraps(func)
            def _func(*args, **kwargs):
                with concurrent.futures.ThreadPoolExecutor(max_workers=1) as \
                        pool:
                    future = pool.submit(func, *args, **kwargs)
                    timed_progress_bar(future,
                                       float_progress,
                                       expected_time,
                                       increments)

                return future.result()

            return _func

        return _progress_bar
