import numpy as np
from westpa.core.states import BasisState, TargetState, InitialState, return_state_type
from westpa.core.segment import Segment
from io import StringIO


class TestStates:
    def test_gen_bstate_from_file(self):
        bstate_txt = """
        # Comment
        0 0.1 a.rst
        1 1e-3 b.rst

        """
        true_answers = [['0', 1e-1, 'a.rst'], ['1', 1e-3, 'b.rst']]

        fake_file = StringIO(bstate_txt)
        bstates = BasisState.states_from_file(fake_file)

        for answer, state in zip(true_answers, bstates):
            assert state.label == answer[0]
            assert state.probability == answer[1]
            assert state.auxref == answer[2]

    def test_gen_tstate_from_file(self):
        tstate_txt = """
        # Comment
        0 1 1
        1 2 3

        """
        true_answers = [['0', np.array([1.0, 1.0])], ['1', np.array([2.0, 3.0])]]

        fake_file = StringIO(tstate_txt)
        tstates = TargetState.states_from_file(fake_file, float)

        for answer, state in zip(true_answers, tstates):
            assert state.label == answer[0]
            assert np.all(state.pcoord == answer[1])

    def test_return_state_type_bstate(self):
        bstate = BasisState(label=0, probability=0.1, state_id=10)
        state_name, state_id = return_state_type(bstate)

        assert state_name == 'BasisState'
        assert state_id == 10

    def test_return_state_type_istate(self):
        istate = InitialState(state_id=2, basis_state_id=10, iter_created=3)
        state_name, state_id = return_state_type(istate)

        assert state_name == 'InitialState'
        assert state_id == 2

    def test_return_state_type_segment(self):
        segment = Segment(seg_id=3)
        state_name, state_id = return_state_type(segment)

        assert state_name == 'Segment'
        assert state_id == 3

    def test_return_state_type_unknown(self):
        state_name, state_id = return_state_type(3)

        assert state_name == 'Unknown'
        assert state_id == float('inf')
