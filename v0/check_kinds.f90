program check_kinds
  print *, "Kind 4: ", kind(0.0_4)
  print *, "Kind 8: ", kind(0.0_8)
  ! print *, "Kind 10: ", kind(0.0_10) ! This might fail if not supported
  ! print *, "Kind 16: ", kind(0.0_16)
end program check_kinds
