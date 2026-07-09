# hocgv-miso: miso-generated Jacobian validity solvers
# License: see hocgv-miso repo

if(TARGET miso)
    return()
endif()

message(STATUS "Third-party: creating target 'miso' via hocgv-miso")

include(CPM)
CPMAddPackage("gh:fsichetti/hocgv-miso#e7c269b4a40156fbdb76db0b0bd20162a3d1c7db")
