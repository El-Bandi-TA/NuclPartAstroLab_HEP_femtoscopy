#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>

std::chrono::time_point<std::chrono::high_resolution_clock> progressbar(
    Int_t event, Int_t nEvents,
    Int_t percentInterval,
    std::chrono::time_point<std::chrono::high_resolution_clock> t_prev,
    Bool_t printBar
) {
    /*Simple progressbar*/

    int currentPercent = event * 100 / nEvents;
    int prevPercent = (event-1) * 100 / nEvents;

    if (!printBar) {
        if (
            event > 1
            && (currentPercent/percentInterval) > (prevPercent/percentInterval)
        ) {
            auto t_now = std::chrono::high_resolution_clock::now();
            auto t_interval = std::chrono::duration_cast<std::chrono::seconds>(
                t_now - t_prev
            ).count();
            t_prev = t_now;

            std::cerr << "Checked " << currentPercent / percentInterval
                      << " interval(s) of " << percentInterval
                      << "% ( " << event << " / " << nEvents << " ) "
                      << "| Last interval took " << t_interval << " s\n"
                      << std::flush;
        }

        return t_prev;
    }

    int barWidth = 70;
    float progress = static_cast<float>(event) / nEvents;

    // print prograss bar
    std::cerr << "\r[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cerr << "=";
        else if (i == pos) std::cerr << ">";
        else std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0)
              << "% ( " << event << " / " << nEvents << " )";
    std::cerr.flush();


    // print timing info
    if (
        event > 1
        && (currentPercent/percentInterval) > (prevPercent/percentInterval)
    ) {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto t_interval = std::chrono::duration_cast<std::chrono::seconds>(
            t_now - t_prev
        ).count();
        t_prev = t_now;

        std::cerr << "\n\r\033[KChecked "
                  << currentPercent / percentInterval << " interval(s) of "
                  << percentInterval << "% | Last interval took " << t_interval
                  << " s\033[A" << std::flush;
    }

    if (currentPercent==100) {std::cerr << "\n\n";}

    return t_prev;
}

std::string get_current_time()
{
    // Get the current time
    auto now = std::chrono::system_clock::now();
    // Convert to time_t for formatting
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    // Convert to tm struct for custom formatting
    std::tm now_tm = *std::localtime(&now_time);
    // Use a string stream to capture formatted output
    std::ostringstream oss;
    oss << std::put_time(&now_tm, "%Y%m%d-%H%M");
    // Return the formatted string
    return oss.str();
}