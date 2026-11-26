#include "RngGoldenAdapter.h"
#include <vector>
#include <type_traits>
#include <cstdio>

using OpenMalaria::LocalRng;
using OpenMalaria::MasterRng;

int main() {
    printf("trivial copy=%d move=%d nothrow move=%d copyable=%d size=%zu\n",
           int(std::is_trivially_copyable_v<LocalRng>), int(std::is_trivially_move_constructible_v<LocalRng>),
           int(std::is_nothrow_move_constructible_v<LocalRng>), int(std::is_copy_constructible_v<LocalRng>),
           sizeof(LocalRng));
#if __has_builtin(__is_trivially_relocatable)
    printf("trivially relocatable=%d\n", int(__is_trivially_relocatable(LocalRng)));
#endif
    for (bool restore : {false, true}) {
        std::vector<LocalRng> rows;
        MasterRng master(0, 0);
        for (int i = 0; i < 4; ++i) {
            rows.emplace_back(master);
            if (restore) {
                std::ostringstream out;
                rows.back().checkpoint(out);
                std::istringstream in(out.str());
                rows.back().checkpoint(in);
            }
            for (auto& rng : rows) {
                gsl_rng view;
                std::memcpy(&view, reinterpret_cast<char*>(&rng) + sizeof(rng) - sizeof(view), sizeof(view));
                printf("restore=%d i=%d engine=%p gsl state=%p type=%p\n", int(restore), i,
                       static_cast<void*>(&rng.m_rng), view.state, static_cast<const void*>(view.type));
                fflush(stdout);
            }
        }
        for (auto& rng : rows) {
            printf("draw=%g\n", rng.gauss(0, 1));
            fflush(stdout);
        }
    }
}
