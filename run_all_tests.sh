#!/bin/bash
(cmake -S . -B test/build && cmake --build test/build) || exit 1 # compile tests
cd test/build

clear

# run tests
if [ $# -eq 0 ]; then
    ./zombit_test && cd ../.. || cd ../..
    exit 0
fi

case $1 in
    input)
        ./zombit_test --gtest_color=yes --gtest_brief=1 --gtest_filter="${2}*" && cd ../.. || cd ../..
        ;;
    *)
        echo ">> nothing tested"
        ;;

esac
