# Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.

add_executable(eepPowerCompare
    eepPowerCompare.cpp
    Settings.cpp
)
target_link_libraries(eepPowerCompare base9 m yaml-cpp)
