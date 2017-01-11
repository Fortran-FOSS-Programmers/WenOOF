
module penf_global_parameters_variables
!-----------------------------------------------------------------------------------------------------------------------------------
!< PENF global (exposed) parameters and variables.
!<
!< @note All module defined entities are public.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
public
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
logical            :: is_initialized = .false. !< Check the initialization of some variables that must be initialized.
integer, parameter :: endianL        = 1       !< Little endian parameter.
integer, parameter :: endianB        = 0       !< Big endian parameter.
integer            :: endian         = endianL !< Bit ordering: Little endian (endianL), or Big endian (endianB).
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Portable kind parameters
#ifdef r16p
integer, parameter :: R16P = selected_real_kind(33,4931) !< 33 digits, range \([10^{-4931}, 10^{+4931} - 1]\); 128 bits.
#else
integer, parameter :: R16P = selected_real_kind(15,307)  !< Defined as R8P; 64 bits.
#endif
integer, parameter :: R8P  = selected_real_kind(15,307)  !< 15 digits, range \([10^{-307} , 10^{+307}  - 1]\); 64 bits.
integer, parameter :: R4P  = selected_real_kind(6,37)    !< 6  digits, range \([10^{-37}  , 10^{+37}   - 1]\); 32 bits.
integer, parameter :: R_P  = R8P                         !< Default real precision.
integer, parameter :: I8P  = selected_int_kind(18)       !< Range \([-2^{63},+2^{63} - 1]\), 19 digits plus sign; 64 bits.
integer, parameter :: I4P  = selected_int_kind(9)        !< Range \([-2^{31},+2^{31} - 1]\), 10 digits plus sign; 32 bits.
integer, parameter :: I2P  = selected_int_kind(4)        !< Range \([-2^{15},+2^{15} - 1]\), 5  digits plus sign; 16 bits.
integer, parameter :: I1P  = selected_int_kind(2)        !< Range \([-2^{7} ,+2^{7}  - 1]\), 3  digits plus sign; 8  bits.
integer, parameter :: I_P  = I4P                         !< Default integer precision.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Format parameters
#ifdef r16p
character(*), parameter :: FR16P  = '(E42.33E4)' !< Output format for kind=R16P real.
#else
character(*), parameter :: FR16P  = '(E23.15E3)' !< Output format for kind=R16P real.
#endif
character(*), parameter :: FR8P   = '(E23.15E3)' !< Output format for kind=R8P real.
character(*), parameter :: FR4P   = '(E13.6E2)'  !< Output format for kind=R4P real.
character(*), parameter :: FR_P   = FR8P         !< Output format for kind=R_P real.
character(*), parameter :: FI8P   = '(I20)'      !< Output format for kind=I8P integer.
character(*), parameter :: FI8PZP = '(I20.19)'   !< Output format for kind=I8P integer with zero prefixing.
character(*), parameter :: FI4P   = '(I11)'      !< Output format for kind=I4P integer.
character(*), parameter :: FI4PZP = '(I11.10)'   !< Output format for kind=I4P integer with zero prefixing.
character(*), parameter :: FI2P   = '(I6)'       !< Output format for kind=I2P integer.
character(*), parameter :: FI2PZP = '(I6.5)'     !< Output format for kind=I2P integer with zero prefixing.
character(*), parameter :: FI1P   = '(I4)'       !< Output format for kind=I1P integer.
character(*), parameter :: FI1PZP = '(I4.3)'     !< Output format for kind=I1P integer with zero prefixing.
character(*), parameter :: FI_P   = FI4P         !< Output format for kind=I_P integer.
character(*), parameter :: FI_PZP = FI4PZP       !< Output format for kind=I_P integer with zero prefixing.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Length (number of digits) of formatted numbers
#ifdef r16p
integer, parameter :: DR16P = 42   !< Number of digits of output format FR16P.
#else
integer, parameter :: DR16P = 23   !< Number of digits of output format FR16P.
#endif
integer, parameter :: DR8P  = 23   !< Number of digits of output format FR8P.
integer, parameter :: DR4P  = 13   !< Number of digits of output format FR4P.
integer, parameter :: DR_P  = DR8P !< Number of digits of output format FR_P.
integer, parameter :: DI8P  = 20   !< Number of digits of output format I8P.
integer, parameter :: DI4P  = 11   !< Number of digits of output format I4P.
integer, parameter :: DI2P  = 6    !< Number of digits of output format I2P.
integer, parameter :: DI1P  = 4    !< Number of digits of output format I1P.
integer, parameter :: DI_P  = DI4P !< Number of digits of output format I_P.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! List of kinds
integer,      parameter :: REAL_KINDS_LIST(1:4)      = [R16P, R8P, R4P, R_P]                    !< List of real kinds.
character(*), parameter :: REAL_FORMATS_LIST(1:4)    = [FR16P, FR8P, FR4P//' ', FR_P]           !< List of real formats.
integer,      parameter :: INTEGER_KINDS_LIST(1:5)   = [I8P, I4P, I2P, I1P,I_P]                 !< List of integer kinds.
character(*), parameter :: INTEGER_FORMATS_LIST(1:5) = [FI8P, FI4P, FI2P//' ', FI1P//' ', FI_P] !< List of integer formats.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Minimum and maximum (representable) values
real(R16P),   parameter :: MinR16P = -huge(1._R16P) !< Minimum value of kind=R16P real.
real(R16P),   parameter :: MaxR16P =  huge(1._R16P) !< Maximum value of kind=R16P real.
real(R8P),    parameter :: MinR8P  = -huge(1._R8P ) !< Minimum value of kind=R8P real.
real(R8P),    parameter :: MaxR8P  =  huge(1._R8P ) !< Maximum value of kind=R8P real.
real(R4P),    parameter :: MinR4P  = -huge(1._R4P ) !< Minimum value of kind=R4P real.
real(R4P),    parameter :: MaxR4P  =  huge(1._R4P ) !< Maximum value of kind=R4P real.
real(R_P),    parameter :: MinR_P  = MinR8P         !< Minimum value of kind=R_P real.
real(R_P),    parameter :: MaxR_P  = MaxR8P         !< Maximum value of kind=R_P real.
integer(I8P), parameter :: MinI8P  = -huge(1_I8P)   !< Minimum value of kind=I8P integer.
integer(I4P), parameter :: MinI4P  = -huge(1_I4P)   !< Minimum value of kind=I4P integer.
integer(I2P), parameter :: MinI2P  = -huge(1_I2P)   !< Minimum value of kind=I2P integer.
integer(I1P), parameter :: MinI1P  = -huge(1_I1P)   !< Minimum value of kind=I1P integer.
integer(I_P), parameter :: MinI_P  = MinI4P         !< Minimum value of kind=I_P integer.
integer(I8P), parameter :: MaxI8P  =  huge(1_I8P)   !< Maximum value of kind=I8P integer.
integer(I4P), parameter :: MaxI4P  =  huge(1_I4P)   !< Maximum value of kind=I4P integer.
integer(I2P), parameter :: MaxI2P  =  huge(1_I2P)   !< Maximum value of kind=I2P integer.
integer(I1P), parameter :: MaxI1P  =  huge(1_I1P)   !< Maximum value of kind=I1P integer.
integer(I_P), parameter :: MaxI_P  =  MaxI4P        !< Maximum value of kind=I_P integer.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Real smallest (representable) values
real(R16P), parameter :: smallR16P = tiny(1._R16P) !< Smallest representable value of kind=R16P real.
real(R8P),  parameter :: smallR8P  = tiny(1._R8P ) !< Smallest representable value of kind=R8P real.
real(R4P),  parameter :: smallR4P  = tiny(1._R4P ) !< Smallest representable value of kind=R4P real.
real(R_P),  parameter :: smallR_P  = smallR8P      !< Smallest representable value of kind=R_P real.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Smallest real representable difference by the running calculator
real(R16P), parameter :: ZeroR16 = nearest(1._R16P, 1._R16P) - &
                                   nearest(1._R16P,-1._R16P) !< Smallest representable difference of kind=R16P real.
real(R8P),  parameter :: ZeroR8  = nearest(1._R8P, 1._R8P) - &
                                   nearest(1._R8P,-1._R8P)   !< Smallest representable difference of kind=R8P real.
real(R4P),  parameter :: ZeroR4  = nearest(1._R4P, 1._R4P) - &
                                   nearest(1._R4P,-1._R4P)   !< Smallest representable difference of kind=R4P real.
real(R_P),  parameter :: Zero    = ZeroR8                    !< Smallest representable difference of kind=R_P real.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Bits/bytes memory requirements (real variables must be computed at runtime)
integer(I2P)            :: BIR16P                         !< Number of bits of kind=R16P real.
integer(I1P)            :: BIR8P                          !< Number of bits of kind=R8P real.
integer(I1P)            :: BIR4P                          !< Number of bits of kind=R4P real.
integer(I1P)            :: BIR_P                          !< Number of bits of kind=R_P real.
integer(I2P)            :: BYR16P                         !< Number of bytes of kind=R16P real.
integer(I1P)            :: BYR8P                          !< Number of bytes of kind=R8P real.
integer(I1P)            :: BYR4P                          !< Number of bytes of kind=R4P real.
integer(I1P)            :: BYR_P                          !< Number of bytes of kind=R_P real.
integer(I8P), parameter :: BII8P = bit_size(MaxI8P)       !< Number of bits of kind=I8P integer.
integer(I4P), parameter :: BII4P = bit_size(MaxI4P)       !< Number of bits of kind=I4P integer.
integer(I2P), parameter :: BII2P = bit_size(MaxI2P)       !< Number of bits of kind=I2P integer.
integer(I1P), parameter :: BII1P = bit_size(MaxI1P)       !< Number of bits of kind=I1P integer.
integer(I_P), parameter :: BII_P = bit_size(MaxI_P)       !< Number of bits of kind=I_P integer.
integer(I8P), parameter :: BYI8P = bit_size(MaxI8P)/8_I8P !< Number of bytes of kind=I8P integer.
integer(I4P), parameter :: BYI4P = bit_size(MaxI4P)/8_I4P !< Number of bytes of kind=I4P integer.
integer(I2P), parameter :: BYI2P = bit_size(MaxI2P)/8_I2P !< Number of bytes of kind=I2P integer.
integer(I1P), parameter :: BYI1P = bit_size(MaxI1P)/8_I1P !< Number of bytes of kind=I1P integer.
integer(I_P), parameter :: BYI_P = bit_size(MaxI_P)/8_I_P !< Number of bytes of kind=I_P integer.
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule penf_global_parameters_variables
