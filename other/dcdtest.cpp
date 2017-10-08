
#ifdef DCDTEST
    // DCD-TEST
    int frame_position;
    frame_position = 0;
    int start,stop,step;

    // start,stop,step
    start = atoi(argv[3]);
    stop = atoi(argv[4]);
    step = atoi(argv[5]);

    /* ---------------------------------------------------------
       Step 3.2 DCD Read.
       --------------------------------------------------------- */
    molfile_timestep_t timestep;
    void *v;
    dcdhandle *dcd;
    int natoms; // from the opening the dcd.
    float sizeMB =0.0, totalMB = 0.0;
    double starttime, endtime, totaltime = 0.0;
    printf("----->  READING DCD  <-----\n");

    // 2. to read a dcd.
    natoms = 0;
    v = open_dcd_read(argv[2],"dcd",&natoms);
    if (!v)
    {
        fprintf(stderr, "main) open_dcd_read failed for file %s\n", *argv);
        return 1;
    }
    dcd = (dcdhandle *)v;
    sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);

    /* ---------------------------------------------------------
       Step 3.3 DCD Load.
       --------------------------------------------------------- */
    frame_position = 1;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1st advance. 1-vmd

    // THIS ONE
    // load_dcd_to_chain(dcd,aa_zero,num_chains);
    // load_dcd_to_atoms(dcd,aa_zero);
    allatoms_0 = load_dcd_to_atoms(dcd,allatoms_0);

    // for(int i=0; i<allatoms_ref.size(); i++)
    // {
    //     std::cout << "xyz-ref: "
    //               << allatoms_ref[i].x << " "
    //               << allatoms_ref[i].y << " "
    //               << allatoms_ref[i].z
    //               << std::endl;
    //     std::cout << "xyz-zer: "
    //               << allatoms_0[i].x << " "
    //               << allatoms_0[i].y << " "
    //               << allatoms_0[i].z
    //               << std::endl;
    // }
    // exit(0);

    frame_position = 2;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 2nd. 2-vmd

    // THIS ONE
    // load_dcd_to_atoms(dcd,aa_later);
    allatoms = load_dcd_to_atoms(dcd,allatoms);
    printf("frame_position: %d\n",frame_position);
    // Advancing Rules.
    // ----------------
    // example. step size -> 5.
    // int advance_size = atoi(argv[3]) - 1;
    // step.
    int advance_size = step - 1; // 0 counts, so advance_size of 4, advances by 5.

    // stop.
    if(stop > dcd->nsets){
        stop = dcd->nsets;
        printf("use stop value: %d\n",stop);
    }

    for (int nset1=2; nset1<start; nset1 += 1 ) {
        // for (int nset1=2; nset1<dcd->nsets; nset1 += step + 1) {
        // debug("forwarding --> current: %d\n",nset1);
        frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);


        // THIS ONE
        // load_dcd_to_chain(dcd,chain_later,num_chains);
        // load_dcd_to_atoms(dcd,aa_later);
        allatoms = load_dcd_to_atoms(dcd,allatoms);

        debug("forwarding --> frame_position: %d\n",frame_position);
    }

    // Get initial starting point.
    printf("--> fast forwarded. to frame: %d\n",frame_position);


    /* ---------------------------------------------------------
       Step 3.4 major for loop begin || doloop.
       --------------------------------------------------------- */
    int nset2;
    nset2 = frame_position;
    do {



        /* ---------------------------------------------------------
           Step 3.6 DCD READ
           --------------------------------------------------------- */
        debug("current: %d\n",nset2);
        printf("frame: --> %d <-- was evaluated.\n",frame_position);

        if (nset2 + advance_size + 1 <= stop)
        {
            frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
            printf("frame: --> %d <-- loaded.\n",frame_position);

            // THIS ONE
            // load_dcd_to_chain(dcd,chain_later,num_chains);
            // load_dcd_to_atoms(dcd,aa_later);
            allatoms = load_dcd_to_atoms(dcd,allatoms);

            // nset2 += advance_size + 1;
        }
        nset2 += advance_size + 1;
    } while (nset2<=stop);

    debug("\n..closing dcd..\n");
    close_file_read(v);
    printf("\n----->  READING DCD completed!  <-----\n");
    printf("\t\tThe maximum possible frame_position was: %d\n",stop);
    printf("\t\tThe last frame evaluated was: %d\n",frame_position);
    exit(0);
#endif // DCDTEST
    // exit(0);
