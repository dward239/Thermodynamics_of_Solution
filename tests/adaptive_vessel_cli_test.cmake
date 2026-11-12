execute_process(COMMAND "${CLI}" "${CASE}" RESULT_VARIABLE result
                OUTPUT_VARIABLE csv ERROR_VARIABLE error)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Adaptive CLI failed: ${error}")
endif()
string(STRIP "${csv}" csv)
string(REPLACE "\n" ";" rows "${csv}")
list(GET rows 0 header)
list(GET rows 1 first)
list(GET rows -1 last)
if(NOT header MATCHES "headspace_H2_mol" OR NOT header MATCHES "solid_mol:XY")
    message(FATAL_ERROR "Missing adaptive or chemistry columns: ${header}")
endif()
string(REPLACE "," ";" initial "${first}")
string(REPLACE "," ";" final "${last}")
list(GET initial 16 initial_daughter)
list(GET final 0 end_time)
list(GET final 3 temperature)
list(GET final 7 headspace_hydrogen)
list(GET final 16 final_daughter)
list(GET final 17 scaled_closure)
list(GET final -1 solid)
if(NOT initial_daughter EQUAL 0 OR NOT end_time EQUAL 10 OR
   temperature LESS 313.205 OR temperature GREATER 313.206 OR
   NOT headspace_hydrogen GREATER 0 OR
   final_daughter LESS 3.011070379e23 OR final_daughter GREATER 3.011070381e23 OR
   scaled_closure GREATER 1 OR NOT solid GREATER 0)
    message(FATAL_ERROR "Unexpected adaptive vessel results: ${first}\n${last}")
endif()

file(READ "${CASE}" model)
string(REGEX REPLACE "chem_component X(.|\n)*$" "" vessel_only "${model}")
file(WRITE "${TEST_DIR}/adaptive_cli_vessel_only.txt" "${vessel_only}")
execute_process(COMMAND "${CLI}" "${TEST_DIR}/adaptive_cli_vessel_only.txt"
                RESULT_VARIABLE result OUTPUT_VARIABLE csv ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT csv MATCHES "headspace_H2_mol" OR csv MATCHES "species_mol:")
    message(FATAL_ERROR "Vessel-only adaptive case failed: ${error}")
endif()

string(REPLACE "323.15" "300" invalid_model "${model}")
set(failed_input "${TEST_DIR}/adaptive_cli_invalid.txt")
set(protected_output "${TEST_DIR}/adaptive_cli_protected.csv")
file(WRITE "${failed_input}" "${invalid_model}")
file(WRITE "${protected_output}" "existing output\n")
execute_process(COMMAND "${CLI}" "${failed_input}" "${protected_output}"
                RESULT_VARIABLE result OUTPUT_VARIABLE csv ERROR_VARIABLE error)
file(READ "${protected_output}" remaining)
if(NOT result EQUAL 1 OR NOT error MATCHES "temperature outside" OR
   NOT remaining STREQUAL "existing output\n")
    message(FATAL_ERROR "An out-of-range case should fail before touching its CSV: ${error}")
endif()
