# hocgv-miso: miso-generated Jacobian validity solvers
# License: see hocgv-miso repo

if(TARGET miso)
    return()
endif()

message(STATUS "Third-party: creating target 'miso' via hocgv-miso")

include(CPM)
CPMAddPackage("gh:fsichetti/hocgv-miso#33973f5bb88db43381c7e11608d37a423e1b105a")
