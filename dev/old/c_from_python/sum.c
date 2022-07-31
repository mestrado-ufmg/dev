/// cc -fPIC -shared -o libsum.so sum.c

// int our_function(int num_numbers) {
//     int i;
//     int sum = 0;
//     for (i = 0; i < num_numbers; i++) {
//         sum += i;
//     }
//     return sum;
// }


float our_function(int num_numbers) {
    int i;
    float sum = 0.1;
    for (i = 0; i < num_numbers; i++) {
        sum += i;
    }
    return sum;
}