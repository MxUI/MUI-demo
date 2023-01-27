#!/bin/bash

DOMAIN_NAME_1="domain1"
DOMAIN_NAME_2="domain2"
DOMAIN_NAME_3="domain3"
INTERFACE_NAME_1="ifs1"
INTERFACE_NAME_2="ifs2"
INTERFACE_NAME_3="ifs3"

mpirun -np 1 multidomain ${DOMAIN_NAME_1} ${INTERFACE_NAME_1} ${INTERFACE_NAME_2} : \
       -np 1 multidomain ${DOMAIN_NAME_2} ${INTERFACE_NAME_2} ${INTERFACE_NAME_3} : \
       -np 1 multidomain ${DOMAIN_NAME_3} ${INTERFACE_NAME_1} ${INTERFACE_NAME_3}
