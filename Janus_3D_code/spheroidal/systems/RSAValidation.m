function RSAValidation()
cubeLength = 10;
Ar = 2;
Deq = 2;
total = 60;

spheroids = RSAAdaptive(cubeLength, Ar, Deq, total);
plotSpheroids(spheroids);

end