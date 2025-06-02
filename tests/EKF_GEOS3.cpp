int main() {
	
	
    eop19620101(21413);
    GGM03S(181);
    DE430Coeff(2285, 1020);
	AuxParamInitialize();
	
	

		
		
	
    int result = all_tests();
	
		
	    std::cout << result << std::endl;
		
		
	if (result == 0)
        printf("PASSED\n");

	
    printf("Tests run: %d\n", tests_run);

    return result != 0;
}