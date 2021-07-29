# Children_ImmersiveVR
Dataset and preprocessing routines accompanying the paper " Immersive virtual reality interferes with default head-trunk coordination strategies in young children".
The study consists of two distinct experiments, which each have their processing pipeline. 

###### Flight Game
The routines Detect_discontinuities.m and createDataTable.m are used consecutively to process the data from this type of recordings
* _Detect_disctontinuities.m_ identifies portions of the continuous kinematic data where "jumps" occured due to some instabilities with the IMU. Portions of the data to exclude are identified as pairs of indices (samples) between which the data is to be excluded and are saved in a text file with similar naming structure as the other files of a recording in the form Sxxx_phase_control_feedback_idxToRemove_timestamp.txt. These files are loaded in the next step. 
* _createDataTable.m_ imports the files containing the kinematic data and flight-related information (trajectory and position of the waypoints), and computes kinematic parameters and performance metrics. The routine saves two data tables, one containing the computed values for all segments, the other one containint average values for each session.

###### Joint Angle Reproduction (JAR) task
* _createDataTable_JAR.m_ imports the files containing the kinematic data and the target positions, and computes kinematic parameters and performance metrics. The routine saves one data table containing average values for each session. 
