caches = None
# ds  = 0.1


size_x = 10
size_y = 10
size_z = 10

RUN_FREE = False
RUN_HALFSPACE = True

volume = (size_y+1) * (size_x+1) * (size_z+1)

RESULTS_FOLDER_FREE = 'results' #'results_box'+str(size_x) +'_'+ str(size_y) + '_'+ str(size_z)
RESULTS_FOLDER_HALFSPACE = 'results_halfspace'
RESULTS_FOLDER_BOX = "results_box_"+str(size_x) +"_"+ str(size_y) + "_"+ str(size_z)


if RUN_FREE:
    RESULTS_FOLDER = RESULTS_FOLDER_FREE
elif  RUN_HALFSPACE:
    RESULTS_FOLDER = RESULTS_FOLDER_HALFSPACE
else:
    RESULTS_FOLDER = RESULTS_FOLDER_BOX