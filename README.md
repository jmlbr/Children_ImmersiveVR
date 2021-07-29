# Children_ImmersiveVR
Dataset and preprocessing routines accompanying the paper " Immersive virtual reality interferes with default head-trunk coordination strategies in young children". The study consists of two distinct experiments investigating the coordination of young children when using immersive VR, which each have their processing pipeline. The protocols are briefly outlined below, for detailed methods please refer to the manuscript. 

#### Study 1

The participants executed an immersive flight game, in which an eagle is steered along a path with coins to collect, either with the head or the torso. The experiments consists of 5 phases: evaluation _Before, Training1, Training 2, evaluation After, evaluation Day After_. The first 4 phases were performed consecutively, the last one on the subsequent day to assess retention. 

#### Study 2

The participants performed a Joint Angle Reproduction test (JAR) done in VR, with the head and the torso. 3 conditions were tested: feedback (a line represents the segment with which the test is done in addition to the target orientation), still (no feedback line), and forward (no feedback line, and a constant forward speed is stimulated). Afterwards, they executed one round of the flight game, steering with their torso, corresponding to the phase _evaluation Before_ of Study 1.


## Matlab code

#### Flight Game
The routines Detect_discontinuities.m and createDataTable.m are used consecutively to process the data from this type of recordings.
* _Detect_discontinuities.m_ identifies portions of the continuous kinematic data where "jumps" occured due to some instabilities with the IMU. Portions of the data to exclude are identified as pairs of indices (samples) between which the data is to be excluded and are saved in a text file with similar naming structure as the other files of a recording in the form Sxxx_phase_control_feedback_idxToRemove_timestamp.txt. These files are loaded in the next step. 
* _createDataTable.m_ imports the files containing the kinematic data and flight-related information (trajectory and position of the waypoints), and computes kinematic parameters and performance metrics. The routine saves two data tables, one containing the computed values for all segments, the other one containint average values for each session.

#### Joint Angle Reproduction (JAR) task
* _createDataTable_JAR.m_ imports the files containing the kinematic data and the target positions, and computes kinematic parameters and performance metrics. The routine saves one data table containing average values for each session. 

The folder _Utils_ contains helper functions. 

## Data
The dataset contains the (raw) data files, subject demographics and processed (Matlab) data tables. 

The folders contain the following files :

#### Flight Game (Studies 1 and 2)
* Raw data files
  - BodyAngles: tridimensional head and torso rotations at each latency,
  - manoeuvreList: list of coin position (waypoints) and the required movement with respect to the preceeding one (baseline = straight forward),
  - TimeAngleRatePosRot : position of the flying object (eagle) at each latency,
  - waypointDistTime: distance to each coin (waypoint), perpendicular to the theoretical flight path; time from previous coin,
  - IdxToRemove: samples containing sensor "jumps", as identified by _Detect_discontinuities.m_.
The files TimeAngleRatePos contain rows filled with '999', which indicate the passing of one coin (specifically crossing the plane containing the coin and perpendicular to the planned flight trajectory). This allows to split the data into trials. The BodyAngle files do not contain these rows; the indices identified in the corresponding TimeAngleRatePos files can be used to split these files. 
* Subject demographics 
* Preprocessed data tables

#### JAR (Study 2)
* Raw data files
  - Body Angles: as described above
  - Events: Time of the display of a new target orientation and corresponding angle
* Subject demographics
* Preprocessed data table

