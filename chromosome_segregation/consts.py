caches = None
# ds  = 0.1


size_x = 100
size_y = 100
size_z = 100

RUN_FREE = False

volume = (size_y+1) * (size_x+1) * (size_z+1)

RESULTS_FOLDER_FREE = 'results' #'results_box'+str(size_x) +'_'+ str(size_y) + '_'+ str(size_z)
RESULTS_FOLDER_HALFSPACE = 'results_halfspace'

if RUN_FREE:
    RESULTS_FOLDER = RESULTS_FOLDER_FREE
else:
    RESULTS_FOLDER = RESULTS_FOLDER_HALFSPACE