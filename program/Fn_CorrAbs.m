function Csr=Fn_CorrAbs(ref,rcv_window)
Csr = abs(sum(ref.*conj(rcv_window)))/length(ref);
end