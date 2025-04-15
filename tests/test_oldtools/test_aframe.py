import sys

import pytest

from westpa.oldtools.aframe.plotting import PlottingMixin


class TestPlotting:
    def test_PlottingMixin(self, monkeypatch):
        a = PlottingMixin()
        with monkeypatch.context() as m:
            # No matplotlib
            m.setattr(a, 'matplotlib_avail', False)

            with pytest.raises(RuntimeError):
                a.require_matplotlib()

        # With Matplotlib
        assert a.require_matplotlib() == sys.modules['matplotlib']
