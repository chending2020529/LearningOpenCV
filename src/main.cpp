/*
 * @Author: chending chending2@seer-group.com
 * @Date: 2023-03-08 10:06:22
 * @LastEditors: chending chending2@seer-group.com
 * @LastEditTime: 2023-03-13 22:04:32
 * @FilePath: /Calibrate2D/src/main.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <fstream>
#include <stdlib.h>

#include "camera.h"
#include "AgvCalibration.h"

int main()
{
    Cameras cam;
    cv::Mat image;
    cam.cameraUSB(0, image);
    // AgvCalibration agv(0.321, "36h11");

    // /*  grab and save image  */
    // // string IP = "rtsp://admin:ld123456@192.168.200.118:554/cam/realmonitor?channel=1&subtype=0";  //
    // // cv::Size board_size = cv::Size(6, 9);
    // // string image_type = "chessboard";
    // // agv.cameraIP(IP, board_size, image_type);
    
    // /*  convert chessboard picture to image(Mat) vector  */
    // int chessboard_num = 15;
    // string chessboard_image_path = "../resource/chessboard/";
    // string chessboard_image_name = "chess_";
    // std::vector<cv::Mat> image_chesses;
    // agv.imageToVector(chessboard_num, chessboard_image_path, chessboard_image_name, image_chesses);
    // std::cout << "size of image vector:" << image_chesses.size() << std::endl;

    // // /*  internal parameters calibration  */
    // cv::Size board_size = cv::Size(6, 9);
    // cv::Size2f board_quare_size = cv::Size2f(36, 36);
    // cv::Mat camera_matrix = cv::Mat(cv::Size(3, 3), CV_32FC1, cv::Scalar(0,0,0));
    // cv::Mat camera_distortion = cv::Mat(cv::Size(5, 1), CV_32FC1, cv::Scalar(0,0,0));;
    // agv.intriCalib(image_chesses, board_size, board_quare_size, camera_matrix, camera_distortion);
    
    // /*  double convert to float  */
    // // agv.cameraMatrixDoubleToFloat(camera_matrix, camera_distortion);

    // // /*  convert pose to image(Mat) vector  */
    // // int pos_num = 14;
    // // string pose_path = "../resource/apriltag/";
    // // string pose_file_name = "pose.txt";
    // // std::vector<cv::Point3f> objPosition;
    // // agv.posToVector(pose_path, pose_file_name, pos_num, objPosition);
    // // std::cout << "size of pos vector:" << objPosition.size() << std::endl;

    // // /*  convert apriltag picture to image(Mat) vector  */
    // // int apriltag_num = 14;
    // // string apriltag_image_path = "../resource/apriltag/";
    // // string apriltag_image_name = "apriltag";
    // // std::vector<cv::Mat> image_apriltag;
    // // agv.imageToVector(apriltag_num, apriltag_image_path, apriltag_image_name, image_apriltag);
    // // std::cout << "size of image vector:" << image_apriltag.size() << std::endl;

    // // /*  internal parameters calibrmapxation  */
    // // std::vector<Point2f> tag_centers;
    // // cv::Mat R = cv::Mat(cv::Size(3,3), CV_32FC1, cv::Scalar::all(0));
    // // cv::Mat T = cv::Mat(cv::Size(3,1), CV_32FC1, cv::Scalar::all(0));
    // // agv.extriCalib(image_apriltag, objPosition, camera_matrix, camera_distortion, tag_centers, R, T);

    // // /*  save camera paremeters  */
    // // string file_name_write = "../resource/result/cameraParameters.xml";
    // // agv.saveParameters(file_name_write, camera_matrix, camera_distortion, R, T);

    // /*  undistortion  */
    // cv::Mat image_input = cv::imread("../resource/chessboard/chess_8.jpg");
    // cv::Size image_size = cv::Size(image_input.cols, image_input.rows);
    // std::cout << "image size:" << image_size << std::endl;
    // cv::Mat image_undistort;
    // cv::Mat R = cv::Mat::eye(3, 3, CV_32F);
    // cv::Mat mapx = cv::Mat(image_size, CV_32FC1);
    // cv::Mat mapy = cv::Mat(image_size, CV_32FC1);

    // /* 求解去畸变后的相机内参 */
    // cv::Rect valid_rect;
    // std::cout << "camera_matrix: " << camera_matrix << std::endl
    //           << "camera_distortion: " << camera_distortion << std::endl; 
    // cv::Mat new_camera_matrix = getOptimalNewCameraMatrix(camera_matrix, camera_distortion, image_size, 1, image_size, &valid_rect, true);
    // std::cout << "new_camera_matrix:" << new_camera_matrix << std::endl;

    // // /* 获取去畸变后矩形框的左上角和右下角 */
    // // cv::Point tl = valid_rect.tl();
    // // int x_undistorted = tl.x;
    // // int y_undistorted = tl.y;

    // // /* 像素平面 >> 归一化平面 */
    // // double fx = new_camera_matrix.at<double>(0, 0);
    // // double fy = new_camera_matrix.at<double>(1, 1);
    // // double cx = new_camera_matrix.at<double>(0, 2);
    // // double cy = new_camera_matrix.at<double>(1, 2);
    // // int x_undistored_normalization = (x_undistorted - cx) / fx;
    // // int y_undistored_normalization = (y_undistorted - cy) / fy;

    // // double k1 = camera_distortion.at<double>(0, 0);
    // // double k2 = camera_distortion.at<double>(0, 1);
    // // double p1 = camera_distortion.at<double>(0, 2);
    // // double p2 = camera_distortion.at<double>(0, 3);
    // // double k3 = camera_distortion.at<double>(0, 4);

    // // /* 求解畸变点 */
    // // double x_distorted, y_distorted;
    // // double r = sqrt(x_distorted * x_distorted + y_distorted * y_distorted);
    // // x_undistored_normalization = x_distorted * (1 + k1 * r * r + k2 * r * r * r * r) + 2 * p1 * x_distorted * y_distorted + p2 * (r * r + 2 * x_distorted * x_distorted);
    // // y_undistored_normalization = y_distorted * (1 + k1 * r * r + k2 * r * r * r * r) + p1 * (r * r + 2 * y_distorted * y_distorted) + 2 * p2 * x_distorted * y_distorted;





    // initUndistortRectifyMap(camera_matrix, camera_distortion, R, new_camera_matrix, image_size, CV_32FC1, mapx, mapy);
    // std::cout << "mapx size: " << mapx.size() << "    mapy size: "<< mapy.size() << std::endl;
    // remap(image_input, image_undistort, mapx, mapy, cv::INTER_LINEAR);
    // // std::cout << "mapx:" << mapx << std::endl;


    // circle(image_input, cv::Point(400, 400), 5, cv::Scalar(0, 0, 255), 2, 8);
    // // float x = mapx.at<float>(400, 400);
    // // float y = mapy.at<float>(400, 400);
    // // // std::cout << "x:" << x << "y:" << y << std::endl;
    // // circle(image_undistort, cv::Point(y, x), 5, cv::Scalar(0, 0, 255), 5, 8);



    // rectangle(image_undistort, valid_rect, cv::Scalar(0, 0, 255), 2, 8);


    // // cv::Point tl = valid_rect.tl();
    // // cv::Point br = valid_rect.br();
    // // std::vector<cv::Point> distort_point;
    // // distort_point.push_back(tl);
    // // distort_point.push_back(br);
    // // std::cout << "distort_point[0]:" << distort_point[0] << "distort_point[1]:" << distort_point[1] << std::endl;
    // // std::vector<cv::Point> undistort_point;
    // // std::cout << "camera_matrix: " << camera_matrix << std::endl;


    // // cv::undistortPoints(distort_point, undistort_point, camera_matrix, camera_distortion);
    // //     std::cout << "1111" << std::endl;
    // // rectangle(image_input, cv::Rect(undistort_point[0], undistort_point[1]), cv::Scalar(0, 0, 255), 2, 8);
    // // cv::undistort(image_input, image_undistort, camera_matrix, camera_distortion);

    // cv::imshow("image_distort", image_input);
    // cv::imshow("image_undistort", image_undistort);
    // cv::imwrite("../resource/result/image_distort8.jpg", image_input);
    // cv::imwrite("../resource/result/image_undistort8.jpg", image_undistort);
    // cv::waitKey(0);












    // /*  pixel point to world point */
    // // std::vector<cv::Point3f> world_points;
    // // double depth = -1.0;
    // // int index = -1;
    // // std::cout << std::endl;
    // // std::cout << "camera_matrix: \n" << camera_matrix <<std::endl;
    // // std::cout << "camera_distortion: \n" << camera_distortion <<std::endl;
    // // std::cout << "R: \n" << R <<std::endl;
    // // std::cout << "T: \n" << T <<std::endl;
    // // agv.ImgPoints2World(tag_centers, camera_matrix, camera_distortion, R, T, world_points, depth, index);
    // // std::cout << "图像坐标： \n" << tag_centers << std::endl;
    // // std::cout << "原始世界坐标： \n" << objPosition << std::endl;
    // // std::cout << "投影世界坐标： \n" << world_points << std::endl;

    return 0;
}