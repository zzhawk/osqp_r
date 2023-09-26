#include "corecrt_math_defines.h"
#include "nlohmann/json.hpp"

// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"

// eigen
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "lane_optimizer.hpp"



int main()
{   
   LaneOptimizer::P_weight w;
   w.l = 1.0;
   w.dl = 20.0;
   w.ddl = 1.0;
   w.dddl = 1.0;
   w.l_end = 1.0;
   w.dl_end = 1.0;
   w.ddl_end = 1.0;
   w.l_ref = 10.0;

   LaneOptimizer::A_setups a;
   a.border = true;
   a.heading = true;
   a.kappa = false;
   a.jerk = false;
   a.start = true;
   a.ddl2dl = true;
   a.dl2l = true;

   auto getLane = []() -> std::unordered_map <std::string, std::vector<double>> {

      std::unordered_map <std::string, std::vector<double>> lane;

      std::ifstream raw("../raw.csv");
      std::string li;
      std::vector<double>buf;
      std::vector<std::string>name = { "S", "l_ref", "l_upper", "l_lower" };
      for (int i = 0; i < name.size(); ++i) {
         lane.emplace(name[i], buf);
      }

      while (std::getline(raw, li)) {
         std::string num = "";
         int i = 0;
         for (auto w : li) {
            if (w != ',') {
               num += w;
            }
            else {
               lane[name[i++]].push_back(std::stod(num));
               num = "";
            }
         }
         lane[name[i++]].push_back(std::stod(num));
      }

      return lane;
   };

   auto getCfg = [&w, &a]() {
      nlohmann::json b;
      std::ifstream cfg("../cfg.json");
      cfg >> b;

      w.l = b.at("weight").at("l");
      w.dl = b.at("weight").at("dl");
      w.ddl = b.at("weight").at("ddl");
      w.dddl = b.at("weight").at("dddl");
      w.l_end = b.at("weight").at("l_end");
      w.dl_end = b.at("weight").at("dl_end");
      w.ddl_end = b.at("weight").at("ddl_end");
      w.l_ref = b.at("weight").at("l_ref");

      a.border = b.at("constrain").at("border");
      a.heading = b.at("constrain").at("heading");
      a.kappa = b.at("constrain").at("kappa");
      a.jerk = b.at("constrain").at("jerk");
      a.start = b.at("constrain").at("start");
      a.ddl2dl = b.at("constrain").at("ddl2dl");
      a.dl2l = b.at("constrain").at("dl2l");
   };

   auto lane = getLane();
   getCfg();

    LaneOptimizer::OptimizationData data;

    auto construct_lane = [&]()->LaneOptimizer::OptimizationData {

        LaneOptimizer::OptimizationData data;

        double ds_last = 0;
        for (int i = 0; i < lane["S"].size(); ++i) {
           double ds = 0.0;
           if (i < lane["S"].size() - 1) ds = lane["S"][i + 1] - lane["S"][i];
           else ds = ds_last;

            SamplePoint sp(ds, lane["l_ref"][i], lane["l_upper"][i], lane["l_lower"][i], 0.0, 0.0);
            data.sample.push_back(sp);
            ds_last = ds;
        }

        data.acc_kappa = 1000.0;
        data.l_start = lane["l_ref"][0];
        data.l_end = lane["l_ref"].back();
        data.dl_start = 0.0;
        data.dl_end = 0.0;
        data.ddl_start = 0.0;
        data.ddl_end = 0.0;
      
        return data;
    };


    LaneOptimizer optmzer(w, a);

    auto ans = optmzer.optimize(construct_lane());

    //for (auto s : ans.s) {
    //   std::cout << s << std::endl;
    //}

    std::ofstream out("../ans.csv");

    for (int i = 0; i < lane["S"].size(); ++i) {
       out << lane["S"][i] << "," <<
          lane["l_ref"][i] << "," <<
          lane["l_upper"][i] << "," <<
          lane["l_lower"][i] << "," <<
          ans.s[i] << "," <<
          std::endl;
    }
}