set(download "https://sciences.ugent.be/skirtextdat/SKIRT9/Resources/")

set(resources Core BPASS AtomsMolecules)
set(Core_expected_version 6)
set(BPASS_expected_version 1)
set(AtomsMolecules_expected_version 3)

foreach (resource IN LISTS resources)
    set(filename SKIRT9_Resources_${resource}_v${${resource}_expected_version})
    file(DOWNLOAD ${download}${filename}.zip "export/resources/${filename}.zip")
    file(ARCHIVE_EXTRACT INPUT "export/resources/${filename}.zip" DESTINATION "export/resources/")
    file(REMOVE "export/resources/${filename}.zip")
    file(COPY "skirt.exe" DESTINATION "export/")
endforeach()