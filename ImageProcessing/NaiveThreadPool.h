// MultiThreadingSandbox.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <functional>
#include <vector>
#include <thread>
#include <condition_variable>
#include <queue>
#include <sstream>
#include <set>

#include "ImageData.h"

class ThreadPool
{
public:
    explicit ThreadPool(std::size_t numThreads)
    {
        start(numThreads);
    }

    ~ThreadPool()
    {
        stop();
    }

    void enqueue(std::function<void()> task)
    {
        {
            std::unique_lock<std::mutex> lock(Mutex);
            TaskQueue.emplace(std::move(task));
        }
        ConditionVar.notify_one();
    }

    bool busy() 
    {
        bool bIsPoolBusy;
        {
            std::unique_lock<std::mutex> lock(Mutex);
            bIsPoolBusy = !WorkingThreadsSet.empty() || !TaskQueue.empty();
        }
        return bIsPoolBusy;
    }

    void stop() noexcept
    {
        {
            std::unique_lock<std::mutex> lock(Mutex);
            bIsStopping = true;
        }

        ConditionVar.notify_all();

        for (std::thread& thread : WorkerThreads)
        {
            thread.join();
        }
    }

private:
    std::vector<std::thread> WorkerThreads;

    std::condition_variable ConditionVar;

    std::mutex Mutex;
    bool bIsStopping = false;

    std::queue<std::function<void()>> TaskQueue;

    std::set<std::thread::id> WorkingThreadsSet;

    void start(std::size_t numThreads)
    {

        for (std::size_t i = 0; i < numThreads; i++)
        {
            WorkerThreads.emplace_back([this]()
            {
                while (true)
                {   
                    // lock queue and get next job
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(Mutex);

                        WorkingThreadsSet.erase(std::this_thread::get_id());
                        ConditionVar.wait(lock, [this]() { return bIsStopping || !TaskQueue.empty();});

                        if (bIsStopping)
                        {
                            break;
                        }

                        WorkingThreadsSet.insert(std::this_thread::get_id());
                        task = std::move(TaskQueue.front());
                        TaskQueue.pop();
                    }
                    // run job
                    task();
                }
            });
        }
    }


};

//int main()
//{
//    {
//        ThreadPool Pool(31);
//        for (int i = 0; i < 100; i++)
//        {
//            Pool.enqueue([i]()
//                {
//                    std::stringstream msg;
//                    msg << i << "\n";
//                    std::cout << msg.str();
//                });
//        }
//
//        while (Pool.busy())
//        {
//
//        }
//    }
//
//    std::cout << "Hello World!\n";
//
//
//
//    std::cout << "Finished\n";
//}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
