set(FOLDER_NAME "merge")

COLMAP_ADD_SOURCES(
 cost_function.h
 similaritytransform.h similaritytransform.cpp
 bundle_adjustment.cpp  bundle_adjustment.h
 bundleadjustment_controller.cpp bundleadjustment_controller.h
 mergeclusters.cpp mergeclusters.h
 
)

