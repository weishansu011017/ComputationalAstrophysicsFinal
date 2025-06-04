#include <vector>
#include <numeric> 
#include <functional>  
#include <cstdlib>
#include <iostream>
#include "OctTree.hpp"

void OctNode::subdivide(std::vector<OctNode>& nodes_list){
    // Center of the current node
    float xmid = 0.5f * (xmin + xmax);
    float ymid = 0.5f * (ymin + ymax);
    float zmid = 0.5f * (zmin + zmax);

    // Get the baseidx
    int baseidx = int(nodes_list.size());

    // Construct 4-new node in the node_list
    nodes_list.emplace_back(xmid, xmax, ymid, ymax, zmax, zmid);    // NEU(0)
    nodes_list.emplace_back(xmin, xmid, ymid, ymax, zmax, zmid);    // NWU(1)
    nodes_list.emplace_back(xmin, xmid, ymin, ymid, zmax, zmid);    // SWU(2)
    nodes_list.emplace_back(xmid, xmax, ymin, ymid, zmax, zmid);    // SEU(3)
    nodes_list.emplace_back(xmid, xmax, ymid, ymax, zmid, zmin);    // NED(4)
    nodes_list.emplace_back(xmin, xmid, ymid, ymax, zmid, zmin);    // NWD(5)
    nodes_list.emplace_back(xmin, xmid, ymin, ymid, zmid, zmin);    // SWD(6)
    nodes_list.emplace_back(xmid, xmax, ymin, ymid, zmid, zmin);    // SED(7)

    // Update childrens
    children[0] = baseidx + 0;
    children[1] = baseidx + 1;
    children[2] = baseidx + 2;
    children[3] = baseidx + 3;
    children[4] = baseidx + 4;
    children[5] = baseidx + 5;
    children[6] = baseidx + 6;
    children[7] = baseidx + 7;

}

void OctTree::reserve_nodes(std::size_t N){
    std::size_t est_nodes = N * 8;

    nodes_list.reserve(est_nodes);
    x.reserve(N);
    y.reserve(N);
    z.reserve(N);
    m.reserve(N);
    order.reserve(N);
}

void OctTree::update_Mtot_COM(){
    for (auto& n : nodes_list) {
        n.Mtot = 0.0f;
        n.COMx = 0.0f;
        n.COMy = 0.0f;
        n.COMz = 0.0f;
    }

    std::function<void(int)> postorder = [&](int idx){
        OctNode& node = nodes_list[idx];

        
        if (node.isLeaf()) {
            float M = 0.0f, Cx = 0.0f, Cy = 0.0f, Cz = 0.0f;
            for (int k = 0; k < node.pcount; ++k) {
                int pid  = order[node.ParticlesLocateidx + k];
                float mj = m[pid];
                M  += mj;
                Cx += mj * x[pid];
                Cy += mj * y[pid];
                Cz += mj * z[pid];
            }
            node.Mtot = M;
            if (M > 0.0f) {
                node.COMx = Cx / M;
                node.COMy = Cy / M;
                node.COMz = Cz / M;
            }
            return;
        }

        float M = 0.0f, Cx = 0.0f, Cy = 0.0f, Cz = 0.0f;
        for (int q = 0; q < 8; ++q) {
            int cidx = node.children[q];
            if (cidx < 0) continue;         
            postorder(cidx);               
            const OctNode& c = nodes_list[cidx];
            M  += c.Mtot;
            Cx += c.Mtot * c.COMx;
            Cy += c.Mtot * c.COMy;
            Cz += c.Mtot * c.COMz;
        }
        node.Mtot = M;
        if (M > 0.0f) {
            node.COMx = Cx / M;
            node.COMy = Cy / M;
            node.COMz = Cz / M;
        }
    };
    postorder(root_idx);
}

void OctTree::build_tree(const std::vector<float>& xin, const std::vector<float>& yin, const std::vector<float>& zin, const std::vector<float>& mass){
    if (xin.size()!=yin.size() || xin.size()!=zin.size() || xin.size()!=mass.size()){
        std::cerr << "x/y/z/mass size mismatch" << std::endl;
        std::exit(1);
    }
    const std::size_t N = xin.size();
    x = xin;
    y = yin;
    z = zin;
    m = mass;
    order.resize(N);
    std::iota(order.begin(), order.end(), 0);

    float xmin = *std::min_element(x.begin(),x.end());
    float xmax = *std::max_element(x.begin(),x.end());
    float ymin = *std::min_element(y.begin(),y.end());
    float ymax = *std::max_element(y.begin(),y.end());
    float zmin = *std::min_element(z.begin(),z.end());
    float zmax = *std::max_element(z.begin(),z.end());
    float eps  = 1e-5f*std::max(std::max(xmax-xmin, ymax-ymin), zmax-zmin);   // prevent partcles fall outside
    xmin-=eps; xmax+=eps; ymin-=eps; ymax+=eps; zmin-=eps; zmax+=eps;

    
    nodes_list.clear();
    nodes_list.reserve(8 * N);
    nodes_list.emplace_back(xmin,xmax,ymin,ymax,zmin,zmax);
    root_idx = 0;
    OctNode& root = nodes_list[0];
    root.ParticlesLocateidx = 0;
    root.pcount = int(N);

    // Build tree
    std::vector<int> stack{root_idx};
    while(!stack.empty()){
        int ni = stack.back(); stack.pop_back();
        OctNode& node = nodes_list[ni];
        if (node.pcount <= leafNmax) continue; // Reach leaf node: Do not need to seperate again

        // ---------- Seperate 8 Octant ----------
        int start = node.ParticlesLocateidx;
        int count = node.pcount;
        float xm  = 0.5f * (node.xmin + node.xmax);
        float ym  = 0.5f * (node.ymin + node.ymax);
        float zm  = 0.5f * (node.zmin + node.zmax);

        // Target: NEU(0) | NWU(1) | SWU(2) | SEU(3) | NED(4) | NWD(5) | SWD(6) | SED(7)
        // Get the first index of this node in the order array
        auto first = order.begin() + start;
        // Get the last index of this node in the order array (start --- count)
        auto last  = first + count;

        // Partition: Move the element that satisfies the condition to front, the others goes back. This function return the mid, which satisfies that first-(satisfy region)-mid-(non-satisfy region)-last
        // First move z
        auto lowerIt = std::partition(first, last, [&](int id){
            return z[id] >= zm;
        });
        // For z > zm move NEU & NWU to front
        auto nenwuEndIt = std::partition(first, lowerIt, [&](int id){
            return y[id] >= ym;
        });
        // Move NEU | NWU
        auto neuEndIt = std::partition(first, nenwuEndIt, [&](int id){
            return x[id] >= xm;
        });
        // Next, deal with SEU & SWU. Move SW|SE 
        auto swuEndIt = std::partition(nenwuEndIt, lowerIt, [&](int id){
            return x[id] < xm;
        });
         // Same process but for D
        auto nenwdEndIt = std::partition(lowerIt, last, [&](int id){
            return y[id] >= ym;
        });
        // Move NED | NWD
        auto nedEndIt = std::partition(lowerIt, nenwdEndIt, [&](int id){
            return x[id] >= xm;
        });
        // Next, deal with SED & SWD. Move SWD|SED 
        auto swdEndIt = std::partition(nenwdEndIt, last, [&](int id){
            return x[id] < xm;
        });

        // Now, the order array is now been seperate into four quadant: order.begin()-##(NE)##-neEndIt-##(NW)##-lowerIt-##(SW)##-swEndIt-##(SE)##-last
        int qptr[8];
        qptr[0] = int(neuEndIt - order.begin());      // NEU
        qptr[1] = int(nenwuEndIt - order.begin());    // NWU
        qptr[2] = int(swuEndIt - order.begin());      // SWU
        qptr[3] = int(lowerIt - order.begin());       // SEU 
        qptr[4] = int(nedEndIt - order.begin());      // NED
        qptr[5] = int(nenwdEndIt - order.begin());    // NWD
        qptr[6] = int(swdEndIt - order.begin());      // SWD
        qptr[7] = int(last - order.begin());          // SED 


        // Generate Node
        int base = int(nodes_list.size());
        node.children[0]=node.children[1]=node.children[2]=node.children[3]=node.children[4]=node.children[5]=node.children[6]=node.children[7]=-1;

        int childStart = start;
        for (int q = 0; q<8; ++q){
            int childCount = qptr[q]-childStart; // Estimate howmuch particles is inside the child node
            if (childCount>0){
                float cx0, cx1, cy0, cy1, cz0, cz1;
                switch(q){   // NEU,NWU,SWU,SEU,NED,NWD,SWD,SED
                    case 0: cx0=xm; cx1=node.xmax; cy0=ym;  cy1=node.ymax; cz0 = zm; cz1 = node.zmax; break;   
                    case 1: cx0=node.xmin; cx1=xm; cy0=ym;  cy1=node.ymax; cz0 = zm; cz1 = node.zmax; break;
                    case 2: cx0=node.xmin; cx1=xm; cy0=node.ymin; cy1=ym; cz0 = zm; cz1 = node.zmax;  break;
                    case 3: cx0=xm; cx1=node.xmax; cy0=node.ymin; cy1=ym; cz0 = zm; cz1 = node.zmax;  break;
                    case 4: cx0=xm; cx1=node.xmax; cy0=ym;  cy1=node.ymax; cz0 = node.zmin; cz1 = zm; break;   
                    case 5: cx0=node.xmin; cx1=xm; cy0=ym;  cy1=node.ymax; cz0 = node.zmin; cz1 = zm; break;
                    case 6: cx0=node.xmin; cx1=xm; cy0=node.ymin; cy1=ym; cz0 = node.zmin; cz1 = zm;  break;
                    default:cx0=xm; cx1=node.xmax; cy0=node.ymin; cy1=ym; cz0 = node.zmin; cz1 = zm;  break;
                }
                nodes_list.emplace_back(cx0, cx1, cy0, cy1, cz0, cz1); // Generate children (Note: Currently `node` variable has been dereference!!)
                OctNode& c = nodes_list.back(); 
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
