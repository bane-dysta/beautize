if(NOT DEFINED EDA_EXECUTABLE OR NOT DEFINED XYZ_INPUT OR NOT DEFINED GJF_INPUT OR
   NOT DEFINED TEST_WORKDIR)
  message(FATAL_ERROR "missing EDA fragment CLI test configuration")
endif()

file(MAKE_DIRECTORY "${TEST_WORKDIR}")

function(run_success label input output)
  execute_process(
    COMMAND "${EDA_EXECUTABLE}" "${input}" ${ARGN} -q -o "${output}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${label} failed\n${stdout}\n${stderr}")
  endif()
  set("${label}_stdout" "${stdout}" PARENT_SCOPE)
endfunction()

run_success(explicit_inline "${XYZ_INPUT}" "${TEST_WORKDIR}/explicit_inline.extxyz"
  --frag water1=1-3 0 --frag water2=4-6 0)
if(NOT explicit_inline_stdout MATCHES "water1[ ]+3[ ]+0" OR
   NOT explicit_inline_stdout MATCHES "water2[ ]+3[ ]+0")
  message(FATAL_ERROR "named --frag definitions were not reported correctly\n${explicit_inline_stdout}")
endif()

run_success(explicit_list "${XYZ_INPUT}" "${TEST_WORKDIR}/explicit_list.extxyz"
  --frag water1=1-3 --frag water2=4-6 --frag-charges 0,0)
run_success(legacy_ids "${XYZ_INPUT}" "${TEST_WORKDIR}/legacy_ids.extxyz"
  --frag-ids 1,1,1,2,2,2 --frag-charges 0,0)
run_success(plain_gjf "${GJF_INPUT}" "${TEST_WORKDIR}/plain_gjf.extxyz"
  --frag water1=1-3 0 --frag water2=4-6 0)
run_success(signed_charges "${XYZ_INPUT}" "${TEST_WORKDIR}/signed_charges.extxyz"
  --frag cation=1-3 1 --frag anion=4-6 -1)
if(NOT signed_charges_stdout MATCHES "cation[ ]+3[ ]+1" OR
   NOT signed_charges_stdout MATCHES "anion[ ]+3[ ]+-1")
  message(FATAL_ERROR "signed inline --frag charges were not parsed correctly\n${signed_charges_stdout}")
endif()

foreach(output IN ITEMS explicit_list_stdout legacy_ids_stdout plain_gjf_stdout)
  if(NOT "${${output}}" MATCHES "Total NCI:[ ]+-2\\.255824 kcal/mol")
    message(FATAL_ERROR "${output} does not reproduce the reference EDA result\n${${output}}")
  endif()
endforeach()

execute_process(
  COMMAND "${EDA_EXECUTABLE}" "${XYZ_INPUT}"
          --frag first=1-3 0 --frag second=3-6 0
          -q -o "${TEST_WORKDIR}/overlap.extxyz"
  RESULT_VARIABLE overlap_result
  OUTPUT_VARIABLE overlap_stdout
  ERROR_VARIABLE overlap_stderr
)
if(overlap_result EQUAL 0 OR
   NOT "${overlap_stdout}${overlap_stderr}" MATCHES "assigned to more than one --frag")
  message(FATAL_ERROR "overlapping --frag selections were not rejected\n${overlap_stdout}\n${overlap_stderr}")
endif()

execute_process(
  COMMAND "${EDA_EXECUTABLE}" "${XYZ_INPUT}"
          --frag first=1-3 0 --frag second=4-5 0
          -q -o "${TEST_WORKDIR}/missing.extxyz"
  RESULT_VARIABLE missing_result
  OUTPUT_VARIABLE missing_stdout
  ERROR_VARIABLE missing_stderr
)
if(missing_result EQUAL 0 OR
   NOT "${missing_stdout}${missing_stderr}" MATCHES "is not assigned by any --frag")
  message(FATAL_ERROR "incomplete --frag coverage was not rejected\n${missing_stdout}\n${missing_stderr}")
endif()

execute_process(
  COMMAND "${EDA_EXECUTABLE}" "${XYZ_INPUT}"
          --frag first=1-3 0 --frag second=4-6
          -q -o "${TEST_WORKDIR}/partial_charge.extxyz"
  RESULT_VARIABLE charge_result
  OUTPUT_VARIABLE charge_stdout
  ERROR_VARIABLE charge_stderr
)
if(charge_result EQUAL 0 OR
   NOT "${charge_stdout}${charge_stderr}" MATCHES "every --frag must provide an integer charge")
  message(FATAL_ERROR "partial inline fragment charges were not rejected\n${charge_stdout}\n${charge_stderr}")
endif()
