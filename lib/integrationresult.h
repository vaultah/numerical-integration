#ifndef INTEGRATIONRESULT_H
#define INTEGRATIONRESULT_H

#include <map>
#include <vector>
#include <memory>
#include "const.h"

namespace Integration {
    template<class Cube>
    struct Result
    {
        std::vector<std::pair<std::shared_ptr<Cube>, REGION_STATE>> cubes;
        long double sum, error;
        std::shared_ptr<Cube> origin;
    };
}

#endif // INTEGRATIONRESULT_H
