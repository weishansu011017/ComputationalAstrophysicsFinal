#include <vector>
#include <numeric> 
#include <functional>  
#include <cstdlib>
#include <iostream>
#include "QuadTree.hpp"

void QuadNode::subdivide(std::vector<QuadNode>& nodes_list){
    // Center of the current node
    float xmid = 0.5f * (xmin + xmax);
    float ymid = 0.5f * (ymin + ymax);

    // Get the baseidx
    int baseidx = int(nodes_list.size());

    // Construct 4-new node in the node_list
    nodes_list.emplace_back(xmid, xmax, ymid, ymax);    // NE(0)
    nodes_list.emplace_back(xmin, xmid, ymid, ymax);    // NW(1)
    nodes_list.emplace_back(xmin, xmid, ymin, ymid);    // SW(2)
    nodes_list.emplace_back(xmid, xmax, ymin, ymid);    // SE(3)

    // Update childrens
    children[0] = baseidx + 0;
    children[1] = baseidx + 1;
    children[2] = baseidx + 2;
    children[3] = baseidx + 3;
}

void QuadTree::reserve_nodes(std::size_t N){
    std::size_t est_nodes = N * 4;

    nodes_list.reserve(est_nodes);
    x.reserve(N);
    y.reserve(N);
    m.reserve(N);
    order.reserve(N);
}

void QuadTree::update_Mtot_COM(){
    for (auto& n : nodes_list) {
        n.Mtot = 0.0f;
        n.COMx = 0.0f;
        n.COMy = 0.0f;
    }

    std::function<void(int)> postorder = [&](int idx){
        QuadNode& node = nodes_list[idx];

        
        if (node.isLeaf()) {
            float M = 0.0f, Cx = 0.0f, Cy = 0.0f;
            for (int k = 0; k < node.pcount; ++k) {
                int pid  = order[node.ParticlesLocateidx + k];
                float mj = m[pid];
                M  += mj;
                Cx += mj * x[pid];
                Cy += mj * y[pid];
            }
            node.Mtot = M;
            if (M > 0.0f) {
                node.COMx = Cx / M;
                node.COMy = Cy / M;
            }
            return;
        }

        float M = 0.0f, Cx = 0.0f, Cy = 0.0f;
        for (int q = 0; q < 4; ++q) {
            int cidx = node.children[q];
            if (cidx < 0) continue;         
            postorder(cidx);               
            const QuadNode& c = nodes_list[cidx];
            M  += c.Mtot;
            Cx += c.Mtot * c.COMx;
            Cy += c.Mtot * c.COMy;
        }
        node.Mtot = M;
        if (M > 0.0f) {
            node.COMx = Cx / M;
            node.COMy = Cy / M;
        }
    };
    postorder(root_idx);
}

void QuadTree::build_tree(const std::vector<float>& xin, const std::vector<float>& yin, const std::vector<float>& mass){
    if (xin.size()!=yin.size() || xin.size()!=mass.size()){
        std::cerr << "x/y/mass size mismatch" << std::endl;
        std::exit(1);
    }
    const std::size_t N = xin.size();
    x = xin;
    y = yin;
    m = mass;
    order.resize(N);
    std::iota(order.begin(), order.end(), 0);

    float xmin = *std::min_element(x.begin(),x.end());
    float xmax = *std::max_element(x.begin(),x.end());
    float ymin = *std::min_element(y.begin(),y.end());
    float ymax = *std::max_element(y.begin(),y.end());
    float eps  = 1e-5f*std::max(xmax-xmin, ymax-ymin);   // prevent partcles fall outside
    xmin-=eps; xmax+=eps; ymin-=eps; ymax+=eps;

    
    nodes_list.clear();
    nodes_list.reserve(4 * N);
    nodes_list.emplace_back(xmin,xmax,ymin,ymax);
    root_idx = 0;
    QuadNode& root = nodes_list[0];
    root.ParticlesLocateidx = 0;
    root.pcount = int(N);

    // Build tree
    std::vector<int> stack{root_idx};
    while(!stack.empty()){
        int ni = stack.back(); stack.pop_back();
        QuadNode& node = nodes_list[ni];
        if (node.pcount <= leafNmax) continue; // Reach leaf node: Do not need to seperate again

        // ---------- Seperate 4 quant ----------
        int start = node.ParticlesLocateidx;
        int count = node.pcount;
        float xm  = 0.5f * (node.xmin + node.xmax);
        float ym  = 0.5f * (node.ymin + node.ymax);

        // Target: NE(0) | NW(1) | SW(2) | SE(3)
        // Get the first index of this node in the order array
        auto first = order.begin() + start;
        // Get the last index of this node in the order array (start --- count)
        auto last  = first + count;

        // Partition: Move the element that satisfies the condition to front, the others goes back. This function return the mid, which satisfies that first-(satisfy region)-mid-(non-satisfy region)-last
        // In this cases. Firt move those id that satisfy (y[id] >= mid (first or second quadant)(NW & NE))
        auto lowerIt = std::partition(first, last, [&](int id){
            return y[id] >= ym;      
        });
        // Next, operating the front part. Move(NE) to front
        auto neEndIt = std::partition(first, lowerIt, [&](int id){
            return x[id] >= xm;  
        });

        // Finally, operating the back part (SW & SE). Move (SE) to back(false)
        auto swEndIt = std::partition(lowerIt, last, [&](int id){
            return x[id] < xm;   
        });

        // Now, the order array is now been seperate into four quadant: order.begin()-##(NE)##-neEndIt-##(NW)##-lowerIt-##(SW)##-swEndIt-##(SE)##-last
        int qptr[4];
        qptr[0] = int(neEndIt - order.begin());   // NE
        qptr[1] = int(lowerIt - order.begin());   // NW
        qptr[2] = int(swEndIt - order.begin());   // SW 
        qptr[3] = int(last - order.begin());   // SE 


        // Generate Node
        int base = int(nodes_list.size());
        node.children[0]=node.children[1]=node.children[2]=node.children[3]=-1;

        int childStart = start;
        for (int q = 0; q<4; ++q){
            int childCount = qptr[q]-childStart; // Estimate howmuch particles is inside the child node
            if (childCount>0){
                float cx0, cx1, cy0, cy1;
                switch(q){   // NE,NW,SW,SE
                    case 0: cx0=xm; cx1=node.xmax; cy0=ym;  cy1=node.ymax; break;
                    case 1: cx0=node.xmin; cx1=xm; cy0=ym;  cy1=node.ymax; break;
                    case 2: cx0=node.xmin; cx1=xm; cy0=node.ymin; cy1=ym;  break;
                    default:cx0=xm; cx1=node.xmax; cy0=node.ymin; cy1=ym;
                }
                nodes_list.emplace_back(cx0,cx1,cy0,cy1); // Generate children (Note: Currently `node` variable has been dereference!!)
                QuadNode& c = nodes_list.back(); 
                c.ParticlesLocateidx = childStart;
                c.pcount = childCount;
                nodes_list[ni].children[q] = (int)nodes_list.size() - 1;
                stack.push_back(nodes_list[ni].children[q]);
            }
            childStart = qptr[q];
        }
    }
    update_Mtot_COM();
}
