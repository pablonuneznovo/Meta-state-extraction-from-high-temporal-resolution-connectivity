function nStates = countNumberOfSequences(sequence)
sequenceLength=length(sequence);
nStates=0;
if sequence(1)==1
    nStates=1;
end
for i=1:sequenceLength-1
    if(sequence(i)==0 && sequence(i+1)==1)
        nStates=nStates+1;
    end
end

end

