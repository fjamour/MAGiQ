/* 
 * File:   simple_logger.h
 * Author: fuad
 *
 * Created on November 13, 2017, 11:47 AM
 */

#ifndef SIMPLE_LOGGER_H
#define	SIMPLE_LOGGER_H

#include <string>
#include <cstdio>
#include <ctime>

using namespace std;

/*
 * Receives logs from all functions, and flushes on demand.
 * TODO must be thread-safe
 */
struct simple_logger_t
{
    string ilog; // holds all logged information
        
    void log
    (
        const string& flog,
        bool print = true,
        bool timed = true
    )
    {
        string str = "";
        if(timed) {
            time_t t = time(0);
            str += ctime(&t);
            str += "  ";
        }
        str += (flog + "\n");
        ilog.append(str);
        if(print)
            printf("%s", str.c_str());
    }
    
    void flush()
    {
        printf("%s\n", ilog.c_str());
    }
};

simple_logger_t LOGGER;


#endif	/* SIMPLE_LOGGER_H */

