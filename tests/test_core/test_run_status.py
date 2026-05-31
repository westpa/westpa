import json

from westpa.core.run_status import SCHEMA_VERSION, read_run_status, status_path_for_h5, write_run_status


class Test_Run_Status:
    def test_status_path_for_h5(self, tmp_path):
        h5file = tmp_path / 'west.h5'

        assert status_path_for_h5(str(h5file)) == str(h5file) + '.progress.json'

    def test_write_and_read_status(self, tmp_path):
        h5file = tmp_path / 'west.h5'

        write_run_status(
            str(h5file),
            {
                'run_state': 'running',
                'phase': 'propagating',
                'current_iteration': 3,
            },
        )
        result = read_run_status(str(h5file))

        assert result.error is None
        assert result.missing is False
        assert result.status['schema_version'] == SCHEMA_VERSION
        assert result.status['west_h5file'] == str(h5file)
        assert result.status['run_state'] == 'running'
        assert result.status['phase'] == 'propagating'
        assert result.status['current_iteration'] == 3

