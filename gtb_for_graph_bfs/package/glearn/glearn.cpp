// gspan.cpp
#include "glearn.h"
#include "lossfunc.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <set>

#include <boost/timer/timer.hpp>

using std::map;
using std::vector;
using std::set;

void Glearn::run(){
    
  // [0] preparation
  double nn = (double) gdata.size();
    
  // [1] setup for the first-level nodes
  map<Triplet,Graph_to_VpairListVector> heap;
  for(unsigned int gid = 0; gid < nn; ++gid){
    VpairList cursor;
    Graph& g = gdata[gid];
    for(unsigned int v=0; v<g.size(); ++v){
      for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){
	if (e->labels.x <= e->labels.z){
	  cursor.set(v,e->to,e->id,0);
	  heap[e->labels][gid].push_back(cursor);
	}
      }
    }
  }

  std::cout << "feat scaling: " << scaling << std::endl;
  std::cout << "# of samples: " << nn << std::endl;
    
  // [2] initialize mu
  for(unsigned int gid = 0; gid < nn; ++gid){ param.mu.push_back(0.0); }
    
  // [3] main iteration of t=1, 2, 3, ...
  for(int rnd=0; rnd<maxitr; ++rnd){
    std::cout << "--------round " << rnd;
    std::cout << " (" << param.nonzero_idx.size() << ") ----------" << std::endl;
        
    // compute the loss with the current theta
    double _loss = 0.0;
    for(unsigned int gid = 0; gid < nn; ++gid){
      _loss += L(gdata[gid].value,param.mu[gid]);
    }
    _loss /= nn;

    double _loglik = _loss;
        
    // compute the penalty terms with the current theta
    double _val;
    for(set<int>::iterator it = param.nonzero_idx.begin(); it != param.nonzero_idx.end(); ++it){
      // before
      //_val = param.theta[s];
      _val = param.theta[*it];
      _loss += lm1 * fabs(_val);
      _loss += 0.5 * lm2 * _val * _val;
    }
    std::cout << "obj.val = " << _loss << " ";
    std::cout << "(loglik=" << _loglik << ", penalty=" << _loss-_loglik << ")" << std::endl;
        
    // initialize the data structure
    param.nonzero_d_idx.clear();
    param.d.clear();
    coc.clear();
        
    param.max_d  = -1.0;
        
    // ------ MAIN SEARCH ------
        
    boost::timer::cpu_timer timer;
    //std::cout << "timer start" << std::endl;
        
    // compute param.d[i] for each pattern i
    loaded_dfscode.resize(1);
    for(map<Triplet,Graph_to_VpairListVector>::iterator it = heap.begin(); it != heap.end(); ++it){
      loaded_dfscode[0].labels = it->first;
      loaded_dfscode[0].time.set(0,1);
      search(it->second,0);
    }
    // -------------------------
        
    std::string result = timer.format();
    //std::cout << "timer end" << result << std::endl;
        
    // compute param.d[0]
    double b0_grad = 0.0;
    double b0_hess = 0.0;
    for(unsigned int gid = 0; gid < nn; ++gid){
      double y = gdata[gid].value;
      b0_grad += d_L_mu(y, param.mu[gid]);
      b0_hess += d2_L_mu(y, param.mu[gid]);
    }
    b0_grad /= nn;
    b0_hess /= nn;
    double H0 = min(max(b0_hess,1e-10),1e10);
    //double z0 = b0_grad - H0*param.theta[0];
    //param.d[0] = -z0/H0 - param.theta[0];

    std::cout << "b0_hess=" << b0_hess << std::endl;
    std::cout << "b0_grad=" << b0_grad << std::endl;

    double d0 = -b0_grad/H0;
    param.d[0] = d0;
    param.delta[0] = b0_grad * d0 + gamma * d0 * d0 * H0;
    param.Hd[0] = fabs(H0*param.d[0]);

    std::cout << "d0=" << d0 << std::endl;

    // remaining elements thta T == 0 but theta != 0
    set<int> result2;
    std::set_difference(param.nonzero_idx.begin(), param.nonzero_idx.end(),
    			param.nonzero_d_idx.begin(), param.nonzero_d_idx.end(),
    			std::inserter(result2, result2.begin()));

    for(set<int>::iterator it=result2.begin(); it!=result2.end(); ++it){
      double val = -param.theta[*it];
      double _hess = 0.0;
      double _grad = 0.0;
      for(set<int>::iterator it2=param.instances[*it].begin(); it2!=param.instances[*it].end(); ++it2){
	double y = gdata[*it2].value;
	_grad += d_L_mu(y,param.mu[*it2])*scaling;
	_hess += d2_L_mu(y,param.mu[*it2])*scaling*scaling;
      }
      _grad /= nn;
      _hess /= nn;

      param.d[*it] = val;
      param.delta[*it]  = (_grad + lm2 * param.theta[*it]) * val + gamma * val * val * (_hess+lm2) + lm1 * fabs(param.theta[*it] + val) - lm1 * fabs(param.theta[*it]);
      param.Hd[*it] = fabs((_hess+lm2)*val);

      // nonzero_d_idx
      if(param.d[*it]!=0){
	param.nonzero_d_idx.insert(*it);
      }

      // update max_d
      if (fabs(val) > param.max_d){
	param.max_d = fabs(val);
      }

      std::cout << "remain [" << *it << "]  d = " << param.d[*it] << std::endl;
    }
        
    // descent direction: Gauss-Southwell-r rule
    //if (rnd % 20 == 0 or rnd < 10){
    //  gs_cnst = max(0.05, 0.95*gs_cnst);
    //}

    double max_Hd = -1.0;
    set<int> deleted;
    for(set<int>::iterator it=param.nonzero_d_idx.begin(); it!=param.nonzero_d_idx.end(); ++it){
      if(fabs(param.d[*it]) < gs_cnst * param.max_d){
	param.d[*it] = 0.0;
	deleted.insert(*it);
      }else if(param.Hd[*it] > max_Hd){
	max_Hd = param.Hd[*it];
      }
    }
    std::cout << "#candidate nonzero d: " << param.nonzero_d_idx.size() << std::endl;
    for(set<int>::iterator it=deleted.begin(); it != deleted.end(); ++it){
      param.nonzero_d_idx.erase(*it);
    }
    std::cout << "#selected nonzero d: " << param.nonzero_d_idx.size() << std::endl;
        
    // ------ CONVERGENCE TEST ------
    std::cout << "max_Hd: " << max_Hd << " ?< " << tol << std::endl;
    std::cout << "max_d: " << param.max_d << std::endl;

    double maxval = max(max_Hd,fabs(b0_grad));
    
    std::cout << "maxval ==> " << maxval << std::endl;

    if (maxval < tol){
      //if (max_Hd < tol){
      std::cout << "CONVERGED: itr:"<< rnd << " feat:" << param.nonzero_idx.size();
      std::cout << " vars:" << param.get_num_vars() << " irreg:" << irregular.size() << std::endl;
      break;
    }

    // recompute delta
    double delta  = param.delta[0];
    for(set<int>::iterator it=param.nonzero_d_idx.begin(); it!=param.nonzero_d_idx.end(); ++it){
      std::cout << "d[" << *it << "] = " << param.d[*it] << std::endl;
      delta += param.delta[*it];
    }

    std::cout << "dirderiv = " << delta << std::endl;
        
    // ------ LINE SEARCH ------
    alpha = min(alpha/pow(c,5),1);
    for(int k = 0; ; ++k){
            
      // for loss
      double loss_before = 0.0;
      double loss_after  = 0.0;
      for(unsigned int gid = 0; gid < nn; ++gid){
      	Graph& g = gdata[gid];
      	double mu_before = param.theta[0];
      	double mu_after  = param.theta[0]+alpha*param.d[0];
      	for(set<int>::iterator it=param.coef_idx[gid].begin(); it!=param.coef_idx[gid].end(); ++it){
      	  mu_before += param.theta[*it]*scaling;
      	  mu_after  += (param.theta[*it]+alpha*param.d[*it])*scaling;
      	}
      	param.mu[gid] = mu_after;
      	loss_before += L(g.value,mu_before);
      	loss_after  += L(g.value,mu_after);
      }
      loss_before /= nn;
      loss_after /= nn;
         
      // for penalty
      double val;
      double penal_before = 0.0;
      double penal_after  = 0.0;
      for(int s=1; s <= (int) param.get_num_vars(); ++s){
      	//for(set<int>::iterator it = param.nonzero_idx.begin(); it != param.nonzero_idx.end(); ++it){
      	//int s = *it;
      	// before
      	val = param.theta[s];
      	penal_before += lm1 * fabs(val);
      	penal_before += 0.5 * lm2 * val * val;
      	// after
      	val = param.theta[s]+alpha*param.d[s];
      	penal_after  += lm1 * fabs(val);
      	penal_after  += 0.5 * lm2 * val * val;
      }
         
      // Armijo rule
      double f_new = loss_after+penal_after;
      double f = loss_before + penal_before;
      std::cout << "  " << k << ") before+delta " << f+alpha*sigma*delta << " >? after " << f_new << std::endl;
      //std::cout << "  " << k << "  before       " << f << "  loss_after" << loss_after << " penal_after " << penal_after << std::endl;
      std::cout << "  " << k << "  before       " << f << std::endl;

      if (f_new-f < alpha*sigma*delta){
	std::cout << "loss: " << loss_after+penal_after << std::endl;
	break;
      }
      alpha = alpha * c;
      if (alpha<1e-4){
	std::cerr << "stepsize too small! Check round off error in fnew-f." << std::endl;
	return;
	//exit(-1);
      }
    }
        
    // ------ PARAMETER UPDATE ------
    
    // update theta0
    double tmp;
    tmp = param.theta[0]+alpha*param.d[0];

    std::cout << std::endl;
    if(param.d[0] != 0){
      std::cout << "c *theta(0) " << tmp << std::endl;
    }else{
      std::cout << "c theta(0) " << tmp << std::endl;
    }
    param.theta[0] = tmp;
        
    // update theta
    int mc = 0;
    int id = 0;
    param.nonzero_idx.clear();
    for(int s=1; s <= (int) param.get_num_vars(); ++s){
      tmp = param.theta[s]+alpha*param.d[s];
      if (tmp != 0.0){
	id++;
	if(param.d[s]==0.0){
	  std::cout << id << "  theta(" << s << ") " << tmp;
	}else{
	  std::cout << id << " *theta(" << s << ") " << tmp;
	}
	std::cout << " [size=" << param.codesize[s] << ",freq=" << param.instances[s].size() << "]" << std::endl;
	param.nonzero_idx.insert(s);
      }
      if (param.d[s]!=0.0) mc++;
      param.theta[s] = tmp;
    }
    std::cout << std::endl;

    std::cout << "nonzero_d:" << param.nonzero_d_idx.size() << std::endl;
    std::cout << "nonzero_d2:" << mc << std::endl;

    // print statistic
    std::cout << "gs_cnst: " << gs_cnst << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "nonzeros: " << param.nonzero_idx.size() << std::endl;
    std::cout << "num vars: " << param.get_num_vars() << std::endl;
    std::cout << "irregulars: " << irregular.size() << std::endl;


    //    std::ostringstream ofstr;
    //    ofstr << ofname << "_" << rnd;

    //    std::ofstream fdump(ofstr.str().c_str());
    //    fdump.setf(std::ios::fixed,std::ios::floatfield);
    //    fdump.precision(12);

    //    fdump << param.theta[0] << std::endl;
    //    for(set<int>::iterator it = param.nonzero_idx.begin(); it != param.nonzero_idx.end(); ++it){
    //      fdump << param.theta[*it] << " " << param.retrieve(*it) << std::endl;
    //    }
    //    
    //    fdump.close();
  }
  // end of main iteration [3]
  //std::cerr << std::endl;
    
}



/////////////////////////

void Glearn::search(Graph_to_VpairListVector& g2tracers, int _parent){

  if ((int) loaded_dfscode.size() > maxpat){ return; }

  if (!is_min_dfs()){ return; }

  // assign the unique ID (pid) to the loaded pattern
  std::ostringstream os;
  os << loaded_dfscode;
  string _code = os.str();

  int pid = param.resister(_code);
  param.codesize[pid] = loaded_dfscode.size();
    
  bool prunable = node_evaluate(g2tracers,pid);

  if (prunable){
    param.h_tmp.erase(pid);
  }else{
    // ----------- visit children ----------- from here
    map<Pair,Graph_to_VpairListVector> b_heap;
    map<int,map<Pair,Graph_to_VpairListVector>,greater<int> > f_heap;
    
    int maxtoc = explore_child_nodes(g2tracers,b_heap,f_heap);
    
    EdgeCode  dcode;
    for(map<Pair,Graph_to_VpairListVector>::iterator it = b_heap.begin(); it != b_heap.end(); ++it){
      dcode.labels = Triplet(-1,it->first.b,-1);
      dcode.time.set(maxtoc, it->first.a);
      loaded_dfscode.push_back(dcode);
      search(it->second,pid);
      loaded_dfscode.pop_back();
    }
      
    for(map<int,map<Pair,Graph_to_VpairListVector>,greater<int> >::iterator
	  it = f_heap.begin(); it != f_heap.end(); ++it){
      for(map<Pair,Graph_to_VpairListVector>::iterator
	    it2 = it->second.begin(); it2 != it->second.end(); ++it2){
	dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
	dcode.time.set(it->first,maxtoc+1);
	loaded_dfscode.push_back(dcode);
	search(it2->second,pid);
	loaded_dfscode.pop_back();
      }
    }
    // ----------- visit children ----------- end here
  }

  // post-order operations....
  postorder_op(pid);
}

bool Glearn::node_evaluate(Graph_to_VpairListVector& g2tracers, int pid){
  // compute grad and hess
  double grad = 0.0;
  double hess = 0.0;
  double mk_upper = 0.0;
  double mk_lower = 0.0;

  // instance bool vector for hashing
  set<int> inst;
    
  // for gid containing the loaded pattern
  for(Graph_to_VpairListVector::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
    int gid = x->first;
    Graph& g = gdata[gid];
    inst.insert(gid);
    param.instances[pid].insert(gid);
    // grad/hess computation
    double tmp = d_L_mu(g.value,param.mu[gid])*scaling;
    grad += tmp; // sc
    double tmp2 = d2_L_mu(g.value,param.mu[gid])*scaling*scaling;
    hess += tmp2; // sc
    // for MK bounds
    if(tmp > 0){
      mk_upper += tmp;
    }else{
      mk_lower += tmp;
    }
  }
    
  double nn = (double) gdata.size();
  grad /= nn;
  hess /= nn;
  mk_upper /= nn;
  mk_lower /= nn;
    
  // irregularity
  unordered_set< set<int> >::iterator irr = coc.find(inst);
    
  bool regular = false;
  if(irr == coc.end()){
    regular = true;
    coc.insert(inst);
  }else{
    irregular.insert(pid);
  }
    
  if(regular || !unique){
    //if(true){
    
    // register instances
    for(set<int>::iterator gid=inst.begin(); gid!=inst.end(); ++gid){
      param.coef_idx[*gid].insert(pid);
    }
    
    // compute H and z
    double theta = param.theta[pid];
    double H  = min(max(hess,1e-10),1e10)+lm2;
    double df = grad + lm2 *theta;
    double z  = df - H*theta;
    // z = df - H*theta = grad + lm2 *theta - (hess+lm2)*theta = grad - hess*theta;
        
    // update formula
    //double ret;
    //if (z < -lm1)    { ret = -(z+lm1)/H; } // ret = -(df-H*theta+lm1)/H = -(df+lm1)/H + theta
    //else if(z > lm1) { ret = -(z-lm1)/H; } // ret = -(df-H*theta-lm1)/H = -(df-lm1)/H + theta
    //else             { ret = 0.0; }
    // compute the descent direction
    //double d = ret - theta;

    double d;
    if (z < -lm1)    { d = -(df+lm1)/H; }
    else if(z > lm1) { d = -(df-lm1)/H; }
    else             { d = -theta; }

    // record d and delta
    param.d[pid] = d;
    param.delta[pid]  = df * d + gamma * d * d * H + lm1 * fabs(theta + d) - lm1 * fabs(theta);
    param.Hd[pid] = fabs(H*d);
    //if(d != 0){
      //std::cout << "  " << pid << ") grad=" << df << " hess=" << H << std::endl;
    //}
        
    // save nonzero_d_idx
    if (d != 0){
      param.nonzero_d_idx.insert(pid);
      // indicate nonzero to all ancestors
      for(map<int,set<int> >::iterator it=param.h_tmp.begin(); it!=param.h_tmp.end(); ++it){
	param.h_tmp[it->first].insert(pid);
      }
    }
      
    // update max_d
    if (fabs(d) > param.max_d){
      param.max_d = fabs(d);
    }
        
  }
    
  param.h_tmp[pid] = set<int>();
    
  set<int> result1;
  std::set_union(param.h1[pid].begin(), param.h1[pid].end(),
		 param.h2[pid].begin(), param.h2[pid].end(),
		 std::inserter(result1, result1.begin()));
  std::set_intersection(result1.begin(), result1.end(),
			param.nonzero_idx.begin(), param.nonzero_idx.end(),
			std::inserter(param.h1[pid], param.h1[pid].begin()));
    
  double b_upper = 0;
  double b_lower = 0;
  //  if(param.h1[pid].size()>0){
  //  std::cout << "map h[" << pid << "] -> ";
  //}
  for(set<int>::iterator it=param.h1[pid].begin(); it!=param.h1[pid].end(); ++it){
    //std::cout << *it << " ";
    double hess = 0.0;
    for(set<int>::iterator it2=param.instances[pid].begin(); it2!=param.instances[pid].end(); ++it2){
      double y = gdata[*it2].value;
      hess += d2_L_mu(y,param.mu[*it2])*scaling*scaling;
    }
    hess /= nn;
    double val = -hess*param.theta[*it]; // (lm2 - H) * th, H = hess + lm2
    if(val > b_upper){
      b_upper = val;
    }else if(val < b_lower){
      b_lower = val;
    }
  }
  //if(param.h1[pid].size()>0){
  //  std::cout << std::endl;
  //}
  return (lm1 >= max(mk_upper+b_upper,-mk_lower-b_lower) && prune);
}

void Glearn::postorder_op(int pid){
  param.h2[pid].clear();
  if(param.h_tmp[pid].size()!=0){
    for(set<int>::iterator it=param.h_tmp[pid].begin(); it!=param.h_tmp[pid].end(); ++it){
      param.h2[pid].insert(*it);
    }
    param.h_tmp.erase(pid);
  }
}

