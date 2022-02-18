import params

# functions that operate on swc arrays
	
def set_data_types(swc_array):

	"""
	Changes data type of id,structure id,radius,parent id, to int
	Changes data type of x,y,z to float

	"""

	cols_int = [
		params.SWC_INDICES['id'],
		params.SWC_INDICES['sid'],
		params.SWC_INDICES['radius'],
		params.SWC_INDICES['pid'],
	]

	cols_float = [
		params.SWC_INDICES['x'],
		params.SWC_INDICES['y'],
		params.SWC_INDICES['z'],
	]

	swc_array[:,cols_int] = swc_array[:,cols_int].astype(float).astype(int)
	swc_array[:,cols_float] = swc_array[:,cols_float].astype(float)
	
	return swc_array

def get_soma_row_from_swc_array(swc_array):
	
	"""
	Gets soma row from swc array

	"""

	soma_row = swc_array[swc_array[:,params.SWC_INDICES['pid']]==params.SOMA_PID]

	if len(soma_row) != 1:
		raise Exception('Error: # of somas is ' + str(len(soma_row)) +' and must be 1!')

	return swc_array[swc_array[:,params.SWC_INDICES['pid']]==params.SOMA_PID][0]


def get_child_rows_from_swc_array(swc_array, parent_id):

	"""
	Gets all the child rows of a given row id from swc array

	"""

	return swc_array[swc_array[:,params.SWC_INDICES['pid']] == parent_id]

