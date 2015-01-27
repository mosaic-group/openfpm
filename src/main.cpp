#include <iostream>
#include "Graph/CartesianGraphFactory.hpp"

#define BOOST_DISABLE_ASSERTS


#define BOOST_TEST_MODULE "C++ test module for OpenFPM_pdata project"
#include <boost/test/included/unit_test.hpp>

#include <grid_dist.hpp>
#include "Point_test.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "memory/HeapMemory.hpp"
#include "Space/Shape/Box.hpp"
#include "util.hpp"

#include "hypercube_unit_test.hpp"
#include "CartesianGraphFactory_unit_test.hpp"
#include "metis_util_unit_test.hpp"
#include "dec_optimizer_unit_test.hpp"
