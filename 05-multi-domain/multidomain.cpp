/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 5                *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
******************************************************************************/

/**
 * @file multidomain.cpp
 * @author Y. H. Tang
 * @date 18 October 2016
 * @brief Test establishing uniface instances for multiple domains
 * using the create_uniface helper function.
 *
 * USAGE mpirun -np n1 ./multidomain domain_name1 interface1 interface2 ... : \
 *              -np n2 ./multidomain domain_name2 interface2 interface3 ... : \
 *              ...
 */

#include "mui.h"

int main( int argc, char ** argv ) {
    mui::mpi_split_by_app();

    std::string domain = argv[1];

    std::vector<std::string> interfaces;
    for ( int i = 2; i < argc; i++ ) interfaces.emplace_back( argv[i] );

    // Option 1: Declare MUI objects using specialisms (i.e. 3 = 3 dimensional, d = double)
    auto ifs = mui::create_uniface<mui::config_3d>( domain, interfaces );

    // Option 2: Declare MUI interface and samplers using templates in config.h
    //auto ifs = mui::create_uniface<mui::default_config>( domain, interfaces );
}
