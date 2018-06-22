    !
    !   Flags for parallel execution
    !   to be cleared on each new execution of parpack
    !
    !
    SUBROUTINE pcontext
    USE PCONTEXT_MODULE
    apps_first = .true.
    aitr_first = .true.
    END SUBROUTINE pcontext