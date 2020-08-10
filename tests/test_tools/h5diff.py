import h5py
from numpy import isclose


class H5Diff:
    def __init__(self, ref_file, test_file):

        self.ref_file = h5py.File(ref_file, 'r')
        self.test_file = h5py.File(test_file, 'r')

    # TODO: I don't love that this is hardcoded to compare the test file against the ref file, but
    #   to fix that I'd have to figure out how to pass an argument when I pass this as a callable to visititems()
    def check_in_ref(self, name, obj):

        # print(f"\nDiffing {name}", end='\n')

        test_object = obj
        ref_object = self.ref_file.get(name)
        assert ref_object is not None

        # If it's a dataset, compare the actual contents
        if type(test_object) == h5py._hl.dataset.Dataset:

            # These are the names of the dtypes corresponding to some h5py Reference objects.
            # Reference objects are pointers to other parts of the h5 file
            # I'm not sure if it's meaningful to directly compare them, I think comparing all the "actual" data
            #   is probably sufficient for now.
            # TODO: Come up with a way to check if two references are equivalent
            ref_names = ['auxref', 'group_ref']

            # This strips the Reference object elements out of the list of elements to compare.
            # HACK: This if statement evaluates to False if the test_object is a single value like a float or something, since
            #   .dtype is referencing an attribute of a numpy array.
            if test_object[()].dtype.names is not None:

                non_reference_test_elements = [
                    x for i, x in enumerate(test_object[()][0]) if not test_object[()].dtype.names[i] in ref_names
                ]

                print(ref_object[()])
                non_reference_ref_elements = [
                    x for i, x in enumerate(ref_object[()][0]) if not ref_object[()].dtype.names[i] in ref_names
                ]

            else:
                non_reference_test_elements = test_object[()]
                non_reference_ref_elements = ref_object[()]

            # This can be either a boolean, as expected from an equality check, or an array of booleans, since if both
            #   things being compared are arrays it becomes an element-wise comparison.
            # If the latter is the case, then we need to ensure ALL elements are true.

            # TODO: This may be checking numbers which differ up to floating point accuracy
            comparison = non_reference_test_elements == non_reference_ref_elements

            if non_reference_ref_elements.dtype == 'float64':
                comparison = isclose(non_reference_test_elements, non_reference_ref_elements)

            if type(comparison) == bool:
                assert non_reference_test_elements == non_reference_ref_elements
            else:
                assert comparison.all()

        # If it's a group, do nothing
        # TODO: Is it sufficient to check only the datasets? The groups should just be organizational units
        elif type(test_object) == h5py._hl.group.Group:
            pass

        # Returning None every loop makes visititems() iterate through every element.
        return None

    # This traverses every element in the test file to see if it exists in the reference.
    # TODO: Is this sensitive to testing A vs B as opposed to B vs A? Need to think about it a bit.
    #   Empirically from some simple testing, seems like it's not
    # Just keep in mind that I'm traversing the elements of the test file, and comparing to the correspondingly named
    #   elements in the reference file.
    # This could be insufficient if the test file just doesn't contain something it's supposed to, which would be a case
    #   for flipping the logic.
    # But if it's flipped, then if there's extra junk in the test file, that wouldn't be revealed.
    # This works for the time being, and as an **extremely** janky hack fix, you could do it both ways...
    def check(self):

        self.test_file.visititems(self.check_in_ref)


# TODO: Add some simple argument handler here
# if __name__ == "__main__":
#
#    diff = H5Diff("./assign_ref.h5", "./analysis_tocompare.h5")
#
#    diff.check()
