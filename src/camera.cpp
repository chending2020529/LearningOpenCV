/*
 * @Author: chending2020529 chending529@gmail.com
 * @Date: 2023-03-13 20:56:17
 * @LastEditors: chending2020529 chending529@gmail.com
 * @LastEditTime: 2023-03-13 23:07:04
 * @FilePath: /LearningOpenCV/src/camera.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "camera.h"

inline void getCurrentTime(std::string &current_time)
{
    std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm *ptm = localtime(&time);
    char date[60] = {0};
    sprintf(date, "%d-%02d-%02d-%02d-%02d-%02d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday,
            (int)ptm->tm_hour, (int)ptm->tm_min, (int)ptm->tm_sec);
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

        char key = static_cast<char>(cv::waitKey(10));                                                              //> fps:10
        if (key == 32)                                                                                              //!> 按空格键保存图片
        {
            std::string current_time;
            getCurrentTime(current_time);
            std::string image_name;
            image_name = "../resource/image-" + current_time + ".jpg";
            std::cout << "图片" << ++number << "保存成功......" << std::endl;
            cv::imwrite(image_name, frame);
        }

        if (key == 27)                                                                                              //!> 按esc键退出视频
            break;                                                                                                  //> 关闭视频
    }
    cv::destroyWindow("USB-camera");
}

void Cameras::cameraRTSP(std::string IP, cv::Mat &image)
{
    cv::VideoCapture captrue(IP);
    if (false == captrue.isOpened())
        std::cout << "can not find IP-camera....." << std::endl;
    else
        std::cout << "IP-camera is working....." << std::endl;

    cv::Mat frame;
    int number = 0;
    cv::namedWindow("IP-camera");
    while (true)
    {
        captrue >> frame;
        cv::imshow("IP-camera", frame);

        char key = static_cast<char>(cv::waitKey(10));                                                             // fps:10
        if(key == 32)                                                                                              //!> 按空格键保存图片
        {
            std::string current_time;
            getCurrentTime(current_time);
            std::string image_name;
            image_name = "../resource/image-" + current_time + ".jpg";
            std::cout << "图片保存成功......" << std::endl;
            cv::imwrite(image_name, frame);
        }

        if (key == 27)                                                                                             //!> 按esc键退出视频
            break;                                                                                                 //> 关闭视频
    }

    cv::destroyWindow("IP-camera");
}