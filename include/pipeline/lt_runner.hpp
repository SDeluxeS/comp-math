#pragma once

#include <string>

namespace pipeline {
    int runLightTime(
        const std::string& metaKernel,
        const std::string& eopFile,
        const std::string& stationFile,
        const std::string& spacecraftFile,
        const std::string& txStationId,
        const std::string& rxStationId,
        double t3TdbSeconds,
        int spacecraftCenterNaif,
        double spacecraftPositionScale,
        const std::string& outputFile
    );
}

int compareWithReference(
    const std::string& metaKernel,
    const std::string& eopFile,
    const std::string& stationFile,
    const std::string& spacecraftFile,
    const std::string& referenceFile,
    int spacecraftCenterNaif,
    double spacecraftPositionScale
);