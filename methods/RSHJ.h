#include "qalsh.h"
#include "heap.h"
#include <queue>
#include <fstream>

class RSHJ
{
public:
    RSHJ(int z, float radius, std::string filepath): z(z), radius(radius) {
        std::ifstream fin(filepath, std::ios::binary);
        fin.read((char*)&n, sizeof(int));
        fin.read((char*)&d, sizeof(int));
        dataset = new float[(uint64_t)n*d];
        fin.read((char*)dataset, (uint64_t)n*d*sizeof(float));
        fin.close();
    }
    ~RSHJ() {
        if (lsh) delete lsh;
        delete[] dataset;
    }

    void preprocess() {
        candidate_list.resize(n);
        r_candidate_list.resize(n);
        lsh = new nns::QALSH<float>(n, d, 2, 0, 2, dataset);
        for (int i = 0; i < n; ++i) {
            candidate_list[i] = lsh->probe(&dataset[(uint64_t)i*d]);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < candidate_list[i].size(); ++j) {
                r_candidate_list[candidate_list[i][j]].push_back(i);
            }
        }
    }
    
    void GSC() {
        color.resize(n);
        cov_freq.resize(n);
        neighbors.resize(n);
        updateable_heap<int, int, std::greater<int>> heap(n + 1);
        for (int i = 0; i < n; ++i) {
            heap.add({i, candidate_list[i].size()});
        }
        while (!heap.empty()) {
            auto top = heap.pop();
            int u = top.first;
            if (color[u] > 0) continue;
            color[u] = 2;
            for (int i = 0; i < r_candidate_list[u].size(); ++i) {
                int v = r_candidate_list[u][i];
                if (color[v] == 0) {
                    heap.minus1(v);
                }
            }
            qs.push_back(u);
            for (int i = 0; i < candidate_list[u].size(); ++i) {
                int v = candidate_list[u][i];
                auto dist = nns::calc_l2_sqr<float>(d, radius, &dataset[(uint64_t)u*d], &dataset[(uint64_t)v*d]);
                if (dist > radius) continue;
                neighbors[u].push_back(v);
                cov_freq[v]++;
                if (color[v] == 0) {
                    color[v] = 1;
                    for (int j = 0; j < r_candidate_list[v].size(); ++j) {
                        int w = r_candidate_list[v][j];
                        if (color[w] == 0) {
                            heap.minus1(w);
                        }
                    }
                }
            }
        }
    }

    void BMC() {
        updateable_heap<int, float, std::less<int>> heap(n + 1);
        for (int i = 0; i < n; ++i) {
            if (color[i] == 1 && cov_freq[i] < z) {
                float sd = 0;
                for (int j = 0; j < n; ++j) {
                    sd += SQR(cov_freq[j] - z);
                }
                for (int j = 0; j < candidate_list[i].size(); ++j) {
                    int v = candidate_list[i][j];
                    if (cov_freq[v] < z) {
                        sd += (cov_freq[v] - z) * 2 + 1;               
                    }
                }
                heap.add({i, sd});
            }
        }
        while (!heap.empty()) {
            auto top = heap.pop();
            int u = top.first;
            if (cov_freq[u] >= z) continue;
            qs.push_back(u);
            color[u] = 2;
            for (int i = 0; i < candidate_list[u].size(); ++i) {
                int v = candidate_list[u][i];
                auto dist = nns::calc_l2_sqr<float>(d, radius, &dataset[(uint64_t)u*d], &dataset[(uint64_t)v*d]);
                if (dist > radius) continue;
                neighbors[u].push_back(v);
                cov_freq[v]++;
            }
            for (int i = 0; i < n; ++i) {
                if (color[i] == 1 && cov_freq[i] < z) {
                    float sd = 0;
                    for (int j = 0; j < n; ++j) {
                        sd += SQR(cov_freq[j] - z);
                    }
                    for (int j = 0; j < candidate_list[i].size(); ++j) {
                        int v = candidate_list[i][j];
                        if (cov_freq[v] < z) {
                            sd += (cov_freq[v] - z) * 2 + 1;               
                        }
                    }
                    heap.update({i, sd});
                }
            }
        }
        std::sort(qs.begin(), qs.end());
        for (auto u : qs) {
            std::sort(neighbors[u].begin(), neighbors[u].end());
        }
    }

    void filter() {
        computed.resize(n);
        int count = 0;
        for (int u : qs) {
            if (neighbors[u].size() < 2) continue;
#pragma omp parallel for sechdule(dynamic) reduction(+:count)
            for (int i = 0; i < neighbors[u].size() - 1; ++i) {
                for (int j = i + 1; j < neighbors[u].size(); ++j) {
                    if (std::find(computed[neighbors[u][i]].begin(), computed[neighbors[u][i]].end(), neighbors[u][j]) != computed[neighbors[u][i]].end()) continue;
                    auto dist = nns::calc_l2_sqr<float>(d, radius, &dataset[(uint64_t)neighbors[u][i]*d], &dataset[(uint64_t)neighbors[u][j]*d]);
                    computed[neighbors[u][i]].push_back(neighbors[u][j]);
                    if (dist <= radius) {
                        count += 2;
                    }
                }
            }
        }
        std::cout << count << std::endl;
    }

    int n;
    int d;
    float *dataset;
    nns::QALSH<float> *lsh;
    std::vector<std::vector<int>> candidate_list;
    std::vector<std::vector<int>> r_candidate_list;

    int z;
    float radius;

    std::vector<int> qs;
    std::vector<int> color;
    std::vector<int> cov_freq;
    std::vector<std::vector<int>> neighbors;
    std::vector<std::vector<int>> computed;
};
