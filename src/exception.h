//
// Created by Egor on 8/11/2024.
//

#ifndef SOOT_AFM_EXCEPTION_H
#define SOOT_AFM_EXCEPTION_H

#include <string>
#include <exception>

class DemException : public std::exception {
public:
    explicit DemException(std::string const & message) : message{message} {}

    const char * what() const noexcept override {
        return message.c_str();
    }

private:
    const std::string message;
};

#endif //SOOT_AFM_EXCEPTION_H
