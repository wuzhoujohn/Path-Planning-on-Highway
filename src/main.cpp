#include <fstream>

#include <math.h>

#include <uWS/uWS.h>

#include <chrono>

#include <iostream>

#include <thread>

#include <vector>

#include "Eigen-3.3/Eigen/Core"

#include "Eigen-3.3/Eigen/QR"

#include "json.hpp"

#include "spline.h"

//Problems
// 1. Can suddenly change mind b/w lane change or not and that doesn't happen smoothly
//2. At one point track vanishes
// 3.
using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() {
  return M_PI;
}
double deg2rad(double x) {
  return x * pi() / 180;
}
double rad2deg(double x) {
  return x * 180 / pi();
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double absolute(double x1, double y1) {
  return sqrt((x1) * (x1) + (y1) * (y1));
}

int mylane(double car_d) {
  int lane_no;
  if ((car_d >= 0) and(car_d < 4.0)) {
    lane_no = 1;
  } else if ((car_d >= 4.0) and(car_d < 8.0)) {
    lane_no = 2;
  } else if ((car_d >= 8.0) and(car_d < 12.0)) {
    lane_no = 3;
  } else {
    lane_no = 0; // out of bound
  }

  return lane_no;

}

int ClosestWaypoint(double x, double y, vector < double > maps_x, vector < double > maps_y) {

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector < double > maps_x, vector < double > maps_y) {

  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y - y), (map_x - x));

  double angle = abs(theta - heading);

  if (angle > pi() / 4) {
    closestWaypoint++;
  }

  return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector < double > getFrenet(double x, double y, double theta, vector < double > maps_x, vector < double > maps_y) {
  int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }

  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;

  double frenet_d = distance(x_x, x_y, proj_x, proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {
    frenet_s,
    frenet_d
  };

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector < double > getXY(double s, double d, vector < double > maps_s, vector < double > maps_x, vector < double > maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp + 1) % maps_x.size();

  double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading - pi() / 2;

  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);

  return {
    x,
    y
  };

}

double mph_to_ms(double mph) {
  return mph * (1600.0 / 3600.0);
}

double ms_to_mph(double ms) {
  return ms * (3600.0 / 1600.0);
}

// Cost Functions
double cost_keep_lane(double dist_closest_front, double num_cars_mylane, double buffer_collision, double cost_collision, int current_lane){
  double cost;
  if (num_cars_mylane == 0) // If my lane is empty
  {
    cost = 0;
    }
	else {
		if (dist_closest_front <= buffer_collision) {
			cost = cost_collision; // too close, set cost to collision cost
  }
		else {
			//car should stay in the middle lane, hence assign a slightly lower cost for staying in middle
			if (current_lane == 1) {
				cost = 50;
			}
			else {
				cost = 100;
			}
			if (dist_closest_front >= 160) {
				// If closest car in front is more than 160 meters away, keep cost as-is
				cost += 0;
			}
			else {
				// If closest car in front is less than 160 meters away, increase cost by 100 
				cost += 100;
			}
		}
	}
  return cost;

}

double cost_change_left(double dist_closest_leftfront, double dist_closest_leftback, double buffer_lc, double cost_collision, double cost_left_turn, double violate_left)

{
  double cost;

  if (violate_left == 1) // Left turn not allowed - collision likely
  {
    cost = cost_collision;
  } else {
	  //if its safe to turn left
	  if ((dist_closest_leftfront >= buffer_lc) and (dist_closest_leftback >= buffer_lc)) // No car is close by on either front or back
    {
		  if (dist_closest_leftfront >= 130) {
			  //even safer to turn left
			  cost = cost_left_turn * 0.6;
		  }
		  else {
        cost = cost_left_turn; // Cost is that of turning left
      }
	  }
	  else {
      cost = cost_collision; // Collision likely
    }
  }

  return cost;

}

double cost_change_right(double dist_closest_rightfront, double dist_closest_rightback, double buffer_lc, double cost_collision, double cost_right_turn, double violate_right)

{
  double cost;

  if (violate_right == 1) // Right turn not allowed - collision likely
  {
    cost = cost_collision;
  } else {
    if ((dist_closest_rightfront >= buffer_lc) and(dist_closest_rightback >= buffer_lc)) // No car is close by on either front or back
    {
		if (dist_closest_rightfront >= 130) {
			//even safer to turn left
			cost = cost_right_turn * 0.6;
		}
		else {
        cost = cost_right_turn; // Cost is that of turning left
      }

    } else {
      cost = cost_collision; // Collision likely
    }
  }

  return cost;

}

double cost_break(double cost_slow_down) {
  double cost;
  cost = cost_slow_down + 0.0;
  return cost;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector < double > map_waypoints_x;
  vector < double > map_waypoints_y;
  vector < double > map_waypoints_s;
  vector < double > map_waypoints_dx;
  vector < double > map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream:: in );

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  // correct spline calculation at the end of track
  map_waypoints_x.push_back(map_waypoints_x[0]);
  map_waypoints_y.push_back(map_waypoints_y[0]);
  map_waypoints_s.push_back(max_s + map_waypoints_s[0]);
  map_waypoints_dx.push_back(map_waypoints_dx[0]);
  map_waypoints_dy.push_back(map_waypoints_dy[0]);

  map_waypoints_x.push_back(map_waypoints_x[1]);
  map_waypoints_y.push_back(map_waypoints_y[1]);
  map_waypoints_s.push_back(max_s + map_waypoints_s[1]);
  map_waypoints_dx.push_back(map_waypoints_dx[1]);
  map_waypoints_dy.push_back(map_waypoints_dy[1]);

  //Decided to fit splines relative to the s coordinate
  //Idea from slack channel 

  tk::spline path_spline_x;
  path_spline_x.set_points(map_waypoints_s, map_waypoints_x);

  tk::spline path_spline_y;
  path_spline_y.set_points(map_waypoints_s, map_waypoints_y);

  tk::spline path_spline_dx;
  path_spline_dx.set_points(map_waypoints_s, map_waypoints_dx);

  tk::spline path_spline_dy;
  path_spline_dy.set_points(map_waypoints_s, map_waypoints_dy);

  h.onMessage([ & map_waypoints_x, & map_waypoints_y, & map_waypoints_s, & map_waypoints_dx, & map_waypoints_dy, & max_s, & path_spline_x, & path_spline_y, & path_spline_dx, & path_spline_dy](uWS::WebSocket < uWS::SERVER > ws, char * data, size_t length,
    uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get < string > ();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];
			int path_size = previous_path_x.size();

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          //Start looking at sensor fusion data

          // cout<<"Sensor Fusion data"<< sensor_fusion;
          vector < double > cars_mylane_s;
			//vector < double > cars_myleft_s;
			//vector < double > cars_myleft_d;
			//vector < double > cars_myright_s;
			//vector < double > cars_myright_d;
			//vector < double > cars_mylane_vel;
			//vector < double > cars_myleft_vel;
			//vector < double > cars_myright_vel;

          double target_d;

          double dist_closest_front = 999;
          double closest_front_vel = 0;

          double dist_closest_leftfront = 999;
          double dist_closest_leftback = 999;
			double laternal_dist_closest_leftfront = 999;
          double closest_leftfront_vel = 0;

          double dist_closest_rightfront = 999;
          double dist_closest_rightback = 999;
			double laternal_dist_closest_rightfront = 999;

          double closest_rightfront_vel = 0;
          double violate_left = 0;
          double violate_right = 0;

          //STEP 1: Analyze the senor fusion data and categorize it meaningfully
          // finding closest cars in all 3 lanes, one in front and one in back
          for (int i = 0; i < sensor_fusion.size(); ++i) {
            auto traffic = sensor_fusion[i];
            double traffic_id = traffic[0];
            double traffic_x = traffic[1];
            double traffic_y = traffic[2];
            double traffic_vx = traffic[3];
            double traffic_vy = traffic[4];
            double traffic_s = traffic[5];
            double traffic_d = traffic[6];
            double traffic_vel = sqrt(traffic_vx * traffic_vx + traffic_vy * traffic_vy);
            double traffic_dist = distance(car_x, car_y, traffic_x, traffic_y);

				double dist = traffic_s - car_s;

				//Categorize traffic as being in one of the 3 lanes, and find closest car
            //My lane
            if (mylane(car_d) == mylane(traffic_d)) {
              cars_mylane_s.push_back(traffic_s);
					//car is in front, as traffic_s > car_s
					if (dist > 0) {
						//find minimal distance for cars in front, and get the speed for the vehicle(in m/s)
						if (dist_closest_front > dist) {
							dist_closest_front = dist;
							closest_leftfront_vel = traffic_vel;
						}
					}
            }

            // My left lane
            else if (mylane(traffic_d) > 0 and(mylane(traffic_d) == mylane(car_d) - 1))
            {
					//car is in front
					if (dist > 0) {
						//if distance is less than 15 meters, we record lateral distance
						if (dist < 15) {
							laternal_dist_closest_leftfront = car_d - traffic_d;
						}
						if (dist_closest_leftfront > dist) {
							dist_closest_leftfront = dist;
						}
					}
					else {
						//car is in back, find the cloest car in the back and assign to dist_closest_leftback
						dist_closest_leftback = min(abs(dist), dist_closest_leftback);
					}
            }

            //My Right lane
            else if (mylane(traffic_d) == mylane(car_d) + 1) {
					//car is in front
					if (dist > 0) {
						//if distance is less than 15 meters, we record lateral distance
						if (dist < 15) {
							laternal_dist_closest_rightfront = traffic_d - car_d;
            }
						if (dist_closest_rightfront > dist) {
							dist_closest_rightfront = dist;
          }
              }
					else {
						//car is in back, find the cloest car in the back and assign to dist_closest_leftback
						dist_closest_rightback = min(abs(dist), dist_closest_rightback);
            }
          }

                }

			//cout << "dist_closest_front:" << dist_closest_front << endl;
			if ((car_d - 4.0) < 0)// Left lane change is not possible
          {
            violate_left = 1;
          }

			if ((car_d + 4.0) > 12) // Right lane is not practical
          {
            violate_right = 1;
          }

          //STEP 2: Calculate cost of actions and choose the one with minimum cost
          double buffer_my = 50;
			double buffer_lc = 40;
			double buffer_collision = 30;
          double cost_collision = 500;
			double cost_left_turn = 200;
			double cost_right_turn = 200;
			double cost_slow_down = 250;

          double continue_lane_tc;
          double break_tc;
          double left_turn_tc;
          double right_turn_tc;

          int decision = 999;
          double min_cost = 800;
          vector < double > costs;

			continue_lane_tc = cost_keep_lane(dist_closest_front, cars_mylane_s.size(), buffer_collision, cost_collision, mylane(car_d));
          costs.push_back(continue_lane_tc);

          break_tc = cost_break(cost_slow_down);
          costs.push_back(break_tc);

          left_turn_tc = cost_change_left(dist_closest_leftfront, dist_closest_leftback, buffer_lc, cost_collision, cost_left_turn, violate_left);
          costs.push_back(left_turn_tc);

          right_turn_tc = cost_change_right(dist_closest_rightfront, dist_closest_rightback, buffer_lc, cost_collision, cost_right_turn, violate_right);
          costs.push_back(right_turn_tc);

          cout<<"Cost continue lane" <<continue_lane_tc<<endl;
          cout<<"Cost slow down" <<break_tc<<endl;
          cout<<"Cost left turn" <<left_turn_tc<<endl;
          cout<<"Cost Right turn" <<right_turn_tc<<endl;

          // Find the action with minimum cost

          for (int i = 0; i < costs.size(); ++i)
          {
            double cost1 = costs[i];
            if (cost1 < min_cost) {
              min_cost = cost1;
              decision = i;
            }
          }

          if (decision == 0) {
            cout << "Decision: Stay in lane with max vel" << endl;
          } else if (decision == 1) {
            cout << "Decision: Stay in lane but slow down" << endl;
          } else if (decision == 2) {
            cout << "Decision: Change to left lane" << endl;
          } else {
            cout << "Decision: Change to right lane" << endl;

          }

          // decision = 0 ; stay in current lane at max speed
					// decision = 1; stay in current lane but slow down
					// decision = 2; change to left lane
					// decision = 3 ; change to right lane

					cout << "Min cost: " << min_cost << endl;
					cout << "Decision is:" << decision << endl;
			 cout << "Velocity of closest front: " << closest_front_vel << endl;

					// STEP 3: Once the car has chosen a decision, decide the best s_inc and target_d for that action

					// Mapping of s and velocity
					// 1 mile/h = 0.447 m/s
					// s = 0.40 ; vel = 45 mile/h = 20.1168 m/s; s = 20.1168 * 0.02s = 0.40 meter
					// s = 0.35; vel = 40;
					// s = 0.30; vel = 35
					const double MAX_SPEED = 50; //mph
			 //The jerk threshold for comfort was approximately 0.3-0.9 m/s3
			 const double MAX_ACC = 0.9; // maximum jerk considering passenger comfort
					//Note closest front vel is likely in m/s -- first convert to miles/hr
					if ((decision == 0) or (decision == 999)) {
						//Stay in current lane at max speed
						//car_speed = 50 / 2.24; //convert mile/h to m/s
				if (car_speed < MAX_SPEED / 2.24) {
					//if speed is lower than maximum speed, increase it by max jerk
					car_speed += MAX_ACC;
				}

				// prevent cars in adjacent lanes hit our car through lane changing(car width roughly equals to 3 m)
				if (laternal_dist_closest_leftfront < 3.5 || laternal_dist_closest_rightfront < 3.5) {
					//slowing down
					car_speed -= MAX_ACC;
				}
				
						//determine which lane car should be on
						if (mylane(car_d) == 1) {
							target_d = 2.0;
						}
						else if (mylane(car_d) == 2) {
							target_d = 6.0;
						}
						else if (mylane(car_d) == 3) {
							target_d = 10.0;
						}
						else target_d = 6.0; //if Car is out of bound, move car to middle lane

					} else if (decision == 1) {
						//Stay in current lane but slow down
				car_speed = car_speed - MAX_ACC;
				//cout << "dist_closest_front: " << dist_closest_front << endl;
				//cout << "closest_front_vel: " << closest_front_vel << endl;

				if ((dist_closest_front < buffer_collision) || (car_speed > closest_front_vel)) {
							car_speed = closest_front_vel;
						}
						cout << "car_speed is:" << car_speed << endl;
						if (mylane(car_d) == 1) {
							target_d = 2.0;
						}
						else if (mylane(car_d) == 2) {
							target_d = 6.0;
						}
						else if (mylane(car_d) == 3) {
							target_d = 10.0;
						}
						else target_d = 6.0; //Car out of bound
					} else if (decision == 2) {
						//Change to left lane
						if (mylane(car_d) <= 2) // You can't switch -- another way to maintain sanity
						{
							target_d = 2.0;
						}
						else if (mylane(car_d) == 3) {
							target_d = 6.0;
						}
						else target_d = 6.0; //Car out of bound

					} else if (decision == 3) {
						//Change to right lane

						if (mylane(car_d) == 1) {
							target_d = 6.0;
						}
						else if (mylane(car_d) >= 2) //You can't switch -- another way to maintain sanity
						{
							target_d = 10.0;
						}
						else target_d = 6.0; //Car out of bound

					} else {
						cout << "Something is wrong";
					}

          //STEP 4: Generate Jerk minimizing trajectories for lane change and acceleration/braking
					vector < double > next_x_vals;
					vector < double > next_y_vals;


					vector<double> pre_and_future_x;
					vector<double> pre_and_future_y;

					double pos_x = car_x;
					double pos_y = car_y;
					double pos_s = car_s;
					double pos_d = car_d;
					double yaw_in_rad = deg2rad(car_yaw);
					//check if car speed is greater than limit, convert 50 mile/h to m/s to align with simulator
					
					if (car_speed > MAX_SPEED / 2.24) {
						car_speed = MAX_SPEED / 2.24;
					}
					cout << "car_speed " << car_speed << endl;

					//As suggested in the lecture use previous path
					if (path_size < 2) {
						double prev_pos_x = pos_x - cos(car_yaw);
						double prev_pos_y = pos_y - sin(car_yaw);

						pre_and_future_x.push_back(prev_pos_x);
						pre_and_future_x.push_back(pos_x);

						pre_and_future_y.push_back(prev_pos_y);
						pre_and_future_y.push_back(pos_y);
					}
					else {
						pos_x = previous_path_x[path_size - 1];
						pos_y = previous_path_y[path_size - 1];

						double pos_x2 = previous_path_x[path_size - 2];
						double pos_y2 = previous_path_y[path_size - 2];
						//angle = atan2(pos_y - pos_y2, pos_x - pos_x2);
						yaw_in_rad = atan2(pos_y - pos_y2, pos_x - pos_x2);

						pre_and_future_x.push_back(pos_x2);
						pre_and_future_x.push_back(pos_x);
						
						pre_and_future_y.push_back(pos_y2);
						pre_and_future_y.push_back(pos_y);
						

					}


					// Setting up three target points in the future.
					vector<double> next_wp0 = getXY(car_s + 30, target_d, map_waypoints_s, map_waypoints_x,
						map_waypoints_y);
					vector<double> next_wp1 = getXY(car_s + 60, target_d, map_waypoints_s, map_waypoints_x,
						map_waypoints_y);
					vector<double> next_wp2 = getXY(car_s + 90, target_d, map_waypoints_s, map_waypoints_x,
						map_waypoints_y);


					pre_and_future_x.push_back(next_wp0[0]);
					pre_and_future_x.push_back(next_wp1[0]);
					pre_and_future_x.push_back(next_wp2[0]);

					pre_and_future_y.push_back(next_wp0[1]);
					pre_and_future_y.push_back(next_wp1[1]);
					pre_and_future_y.push_back(next_wp2[1]);


					// Making coordinates to local car coordinates.
					for (int i = 0; i < pre_and_future_x.size(); i++) {
						double shift_x = pre_and_future_x[i] - pos_x; //difference between previous/future x and current x
						double shift_y = pre_and_future_y[i] - pos_y; //difference between previous/future y and current y

						pre_and_future_x[i] = shift_x * cos(0 - yaw_in_rad) - shift_y * sin(0 - yaw_in_rad);
						pre_and_future_y[i] = shift_x * sin(0 - yaw_in_rad) + shift_y * cos(0 - yaw_in_rad);
					}

					// Create the spline.
					tk::spline s;
					s.set_points(pre_and_future_x, pre_and_future_y); // using 2 previous points and 3 furture points to initialize the spline

					// Output path points from previous path for continuity.
					for (int i = 0; i < path_size; i++) {
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					}

					// Calculate distance y position on 30 meters ahead.
					double target_x = 30.0;
					double target_y = s(target_x);
					double target_dist = sqrt(target_x * target_x + target_y * target_y);
					double x_add_on = 0;
					
					//cout << "50 - path_size " << 50 - path_size << endl;
					for (int i = 1; i < 50 - path_size; i++) {
						double N = target_dist / (0.02 * car_speed); // divided into N segments

						double x_point = x_add_on + target_x / N;
						double y_point = s(x_point); // To calculate distance y position .

						x_add_on = x_point;

						double x_ref = x_point;
						double y_ref = y_point;

						//x_point here is the increment of current x coordinate
						x_point = x_ref * cos(yaw_in_rad) - y_ref * sin(yaw_in_rad);
						y_point = x_ref * sin(yaw_in_rad) + y_ref * cos(yaw_in_rad);

						//increce current x by the increment calculated above
						x_point += pos_x;
						y_point += pos_y;

						next_x_vals.push_back(x_point);
						next_y_vals.push_back(y_point);
					}
					
					cout << endl;
					cout << endl;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\"," + msgJson.dump() + "]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse * res, uWS::HttpRequest req, char * data,
    size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res -> end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res -> end(nullptr, 0);
    }
  });

  h.onConnection([ & h](uWS::WebSocket < uWS::SERVER > ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([ & h](uWS::WebSocket < uWS::SERVER > ws, int code,
    char * message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen("127.0.0.1", port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}