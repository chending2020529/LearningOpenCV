/*
 * @Author: chending2020529 chending529@gmail.com
 * @Date: 2023-03-13 20:56:17
 * @LastEditors: chending2020529 chending529@gmail.com
 * @LastEditTime: 2023-03-13 21:58:39
 * @FilePath: /LearningOpenCV/src/camera.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "camera.h"

inline void getCurrentTime(std::string &current_time)
{
    auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    auto local_time_now = std::chrono::system_clock::now().time_since_epoch();
    std::chrono::milliseconds local_time_now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(local_time_now);
    struct tm *ptm = localtime(&time);
    char date[60] = {0};
    sprintf(date, "%d-%02d-%02d-%02d-%02d-%02d-%3d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday,
            (int)ptm->tm_hour, (int)ptm->tm_min, (int)ptm->tm_sec, local_time_now_ms % 1000);
    current_time = std::string(date);
}

void Cameras::cameraUSB(int camera_index, cv::Mat &image)
{
    cv::VideoCapture captrue(camera_index);
    if (false == captrue.isOpened())
        std::cout << "can not find USB-camera....." << std::endl;
    else
        std::cout << "USB-camera is working....." << std::endl;

    cv::Mat frame;
    int number = 0;
    cv::namedWindow("USB-camera");
    while (true)
    {
        captrue >> frame;
        cv::imshow("USB-camera", frame);

        char key = static_cast<char>(cv::waitKey(10)); // fps:10
        // if (key == 32)                                 // click " " to save image
        // {
        //     std::vector<cv::Point2f> image_points;
        //     cv::Size board_size;
        //     bool found_chessboard = cv::findChessboardCorners(frame, board_size, image_points,
        //                                                         cv::CALIB_CB_FAST_CHECK /*| cv::CALIB_CB_ADAPTIVE_THRESH
        //                                                         | cv::CALIB_CB_NORMALIZE_IMAGE*/
        //     );
        //     if (found_chessboard)
        //     {
        //         std::string current_time;
        //         getCurrentTime(current_time);
        //         std::string image_name;
        //         image_name = "../resource/chessboard/chess_" + current_time + std::to_string(++number) + ".jpg";
        //         std::cout << "find corners success.........." << std::endl;
        //         cv::imwrite(image_name, frame);
        //     }
        //     else
        //     {
        //         std::cout << "can not find corners!!!!!! \n";
        //     }
        // }

        if (key == 27) // click "esc" to close video
        {
            break; // 关闭视频
        }
    }
    cv::destroyWindow("video");
}

void Cameras::cameraRTSP(std::string IP, cv::Mat &image)
{
    cv::VideoCapture captrue(IP);
    if (false == captrue.isOpened())
        std::cout << "can not find video....." << std::endl;
    else
        std::cout << "video is working....." << std::endl;

    cv::Mat frame;
    int number = 0;
    cv::namedWindow("video");
    while (true)
    {
        captrue >> frame;
        cv::imshow("video", frame);

        char key = static_cast<char>(cv::waitKey(10)); // fps:10
        // if (key == 32)                                 // click " " to save image
        // {
        //     if (image_type == "chessboard")
        //     {
        //         std::vector<cv::Point2f> image_points;
        //         bool found_chessboard = cv::findChessboardCorners(frame, board_size, image_points,
        //                                                           cv::CALIB_CB_FAST_CHECK);
        //         if (found_chessboard)
        //         {
        //             std::string current_time;
        //             getCurrentTime(current_time);
        //             std::string image_name;
        //             image_name = "../resource/chessboard/chessboard-" + std::to_string(++number) + ".jpg";
        //             std::cout << "picture " << std::to_string(number) << " find corners success.........." << std::endl;
        //             cv::imwrite(image_name, frame);
        //         }
        //         else
        //         {
        //             std::cout << "can not find corners！！！！！ \n";
        //         }
        //     }
        // }

        if (key == 27) // click "esc" to close video
            break;     // 关闭视频
    }

    cv::destroyWindow("video");
}