package OptimalExperimentalDesign;

import java.lang.foreign.Arena;
import java.lang.foreign.AddressLayout;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemoryLayout;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SegmentAllocator;
import java.lang.foreign.SymbolLookup;
import java.lang.foreign.ValueLayout;
import java.lang.foreign.MemoryLayout.PathElement;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

public final class WebGpuPanama {
    private static final Linker LINKER = Linker.nativeLinker();
    private static final Arena LIBRARY_ARENA = Arena.ofShared();
    private static final Arena UPCALL_ARENA = Arena.ofShared();
    private static final Path LIB_PATH = findLibraryPath();
    private static final SymbolLookup LOOKUP = SymbolLookup.libraryLookup(LIB_PATH, LIBRARY_ARENA);

    private static final AddressLayout PTR = ValueLayout.ADDRESS;
    private static final ValueLayout.OfInt I32 = ValueLayout.JAVA_INT;
    private static final ValueLayout.OfLong I64 = ValueLayout.JAVA_LONG;

    private static final long WGPU_STRLEN = -1L;
    private static final int WGPU_TRUE = 1;
    private static final int WGPU_FALSE = 0;
    private static final int WGPUCallbackMode_AllowSpontaneous = 3;
    private static final int WGPUBackendType_Undefined = 0;
    private static final int WGPUPowerPreference_HighPerformance = 2;
    private static final int WGPUFeatureLevel_Undefined = 0;
    private static final int WGPUBufferBindingType_Uniform = 2;
    private static final int WGPUBufferBindingType_Storage = 3;
    private static final int WGPUBufferBindingType_ReadOnlyStorage = 4;
    private static final long WGPUBufferUsage_MapRead = 1L;
    private static final long WGPUBufferUsage_CopySrc = 4L;
    private static final long WGPUBufferUsage_CopyDst = 8L;
    private static final long WGPUBufferUsage_Uniform = 64L;
    private static final long WGPUBufferUsage_Storage = 128L;
    private static final long WGPUMapMode_Read = 1L;
    private static final long WGPUShaderStage_Compute = 4L;
    private static final int WGPURequestAdapterStatus_Success = 1;
    private static final int WGPURequestDeviceStatus_Success = 1;
    private static final int WGPUMapAsyncStatus_Success = 1;
    private static final long WGPU_WHOLE_SIZE = -1L;

    private static final MemoryLayout STRING_VIEW = MemoryLayout.structLayout(
        PTR.withName("data"),
        I64.withName("length")
    ).withByteAlignment(8);
    private static final MemoryLayout SHADER_SOURCE_WGSL = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        I32.withName("sType"),
        MemoryLayout.paddingLayout(4),
        PTR.withName("code_data"),
        I64.withName("code_length")
    );
    private static final MemoryLayout REQUEST_ADAPTER_OPTIONS = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        I32.withName("featureLevel"),
        I32.withName("powerPreference"),
        I32.withName("forceFallbackAdapter"),
        I32.withName("backendType"),
        PTR.withName("compatibleSurface")
    );
    private static final MemoryLayout REQUEST_ADAPTER_CALLBACK_INFO = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        I32.withName("mode"),
        MemoryLayout.paddingLayout(4),
        PTR.withName("callback"),
        PTR.withName("userdata1"),
        PTR.withName("userdata2")
    );
    private static final MemoryLayout REQUEST_DEVICE_CALLBACK_INFO = REQUEST_ADAPTER_CALLBACK_INFO;
    private static final MemoryLayout BUFFER_MAP_CALLBACK_INFO = REQUEST_ADAPTER_CALLBACK_INFO;
    
    private static final MemoryLayout DEVICE_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length"),
        I64.withName("requiredFeatureCount"),
        PTR.withName("requiredFeatures"),
        PTR.withName("requiredLimits"),
        PTR.withName("defaultQueue_nextInChain"),
        PTR.withName("defaultQueue_label_data"),
        I64.withName("defaultQueue_label_length"),
        PTR.withName("deviceLostCallbackInfo_nextInChain"),
        I32.withName("deviceLostCallbackInfo_mode"),
        MemoryLayout.paddingLayout(4),
        PTR.withName("deviceLostCallbackInfo_callback"),
        PTR.withName("deviceLostCallbackInfo_userdata1"),
        PTR.withName("deviceLostCallbackInfo_userdata2"),
        PTR.withName("uncapturedErrorCallbackInfo_nextInChain"),
        PTR.withName("uncapturedErrorCallbackInfo_callback"),
        PTR.withName("uncapturedErrorCallbackInfo_userdata1"),
        PTR.withName("uncapturedErrorCallbackInfo_userdata2")
    );
    private static final MemoryLayout SHADER_MODULE_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length")
    );
    private static final MemoryLayout BUFFER_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length"),
        I64.withName("usage"),
        I64.withName("size"),
        I32.withName("mappedAtCreation"),
        MemoryLayout.paddingLayout(4)
    );
    private static final MemoryLayout PIPELINE_LAYOUT_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length"),
        I64.withName("bindGroupLayoutCount"),
        PTR.withName("bindGroupLayouts"),
        I32.withName("immediateSize")
    );
    private static final MemoryLayout COMPUTE_PIPELINE_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length"),
        PTR.withName("layout"),
        PTR.withName("compute_nextInChain"),
        PTR.withName("compute_module"),
        PTR.withName("compute_entryPoint_data"),
        I64.withName("compute_entryPoint_length"),
        I64.withName("compute_constantCount"),
        PTR.withName("compute_constants")
    );
    private static final MemoryLayout COMPUTE_PASS_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length"),
        PTR.withName("timestampWrites")
    );
    private static final MemoryLayout COMMAND_ENCODER_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        PTR.withName("label_data"),
        I64.withName("label_length")
    );
    private static final MemoryLayout COMMAND_BUFFER_DESCRIPTOR = COMMAND_ENCODER_DESCRIPTOR;
    private static final MemoryLayout WGPU_INSTANCE_DESCRIPTOR = MemoryLayout.structLayout(
        PTR.withName("nextInChain"),
        I64.withName("requiredFeatureCount"),
        PTR.withName("requiredFeatures"),
        PTR.withName("requiredLimits")
    );

    private static final MethodHandle MH_CREATE_INSTANCE = downcall("wgpuCreateInstance",
        FunctionDescriptor.of(PTR, PTR));
    private static final MethodHandle MH_INSTANCE_REQUEST_ADAPTER = downcall("wgpuInstanceRequestAdapter",
        FunctionDescriptor.of(PTR, PTR, PTR, REQUEST_ADAPTER_CALLBACK_INFO));
    private static final MethodHandle MH_ADAPTER_REQUEST_DEVICE = downcall("wgpuAdapterRequestDevice",
        FunctionDescriptor.of(PTR, PTR, PTR, REQUEST_DEVICE_CALLBACK_INFO));
    private static final MethodHandle MH_DEVICE_GET_QUEUE = downcall("wgpuDeviceGetQueue",
        FunctionDescriptor.of(PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_SHADER_MODULE = downcall("wgpuDeviceCreateShaderModule",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_BUFFER = downcall("wgpuDeviceCreateBuffer",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_COMMAND_ENCODER = downcall("wgpuDeviceCreateCommandEncoder",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_BIND_GROUP_LAYOUT = downcall("wgpuDeviceCreateBindGroupLayout",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_PIPELINE_LAYOUT = downcall("wgpuDeviceCreatePipelineLayout",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_COMPUTE_PIPELINE = downcall("wgpuDeviceCreateComputePipeline",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_DEVICE_CREATE_BIND_GROUP = downcall("wgpuDeviceCreateBindGroup",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_QUEUE_WRITE_BUFFER = downcall("wgpuQueueWriteBuffer",
        FunctionDescriptor.ofVoid(PTR, PTR, I64, PTR, I64));
    private static final MethodHandle MH_COMMAND_ENCODER_BEGIN_COMPUTE_PASS = downcall("wgpuCommandEncoderBeginComputePass",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_COMPUTE_PASS_SET_PIPELINE = downcall("wgpuComputePassEncoderSetPipeline",
        FunctionDescriptor.ofVoid(PTR, PTR));
    private static final MethodHandle MH_COMPUTE_PASS_SET_BIND_GROUP = downcall("wgpuComputePassEncoderSetBindGroup",
        FunctionDescriptor.ofVoid(PTR, I32, PTR, I64, PTR));
    private static final MethodHandle MH_COMPUTE_PASS_DISPATCH = downcall("wgpuComputePassEncoderDispatchWorkgroups",
        FunctionDescriptor.ofVoid(PTR, I32, I32, I32));
    private static final MethodHandle MH_COMPUTE_PASS_END = downcall("wgpuComputePassEncoderEnd",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_COMMAND_COPY_BUFFER_TO_BUFFER = downcall("wgpuCommandEncoderCopyBufferToBuffer",
        FunctionDescriptor.ofVoid(PTR, PTR, I64, PTR, I64, I64));
    private static final MethodHandle MH_COMMAND_ENCODER_FINISH = downcall("wgpuCommandEncoderFinish",
        FunctionDescriptor.of(PTR, PTR, PTR));
    private static final MethodHandle MH_QUEUE_SUBMIT_FOR_INDEX = downcall("wgpuQueueSubmitForIndex",
        FunctionDescriptor.of(I64, PTR, I64, PTR));
    private static final MethodHandle MH_BUFFER_MAP_ASYNC = downcall("wgpuBufferMapAsync",
        FunctionDescriptor.of(I64, PTR, I64, I64, I64, BUFFER_MAP_CALLBACK_INFO));
    private static final MethodHandle MH_BUFFER_GET_CONST_MAPPED_RANGE = downcall("wgpuBufferGetConstMappedRange",
        FunctionDescriptor.of(PTR, PTR, I64, I64));
    private static final MethodHandle MH_BUFFER_UNMAP = downcall("wgpuBufferUnmap",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_DEVICE_POLL = downcall("wgpuDevicePoll",
        FunctionDescriptor.of(I32, PTR, I32, PTR));
    private static final MethodHandle MH_RELEASE = downcall("wgpuBufferRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_BIND_GROUP_LAYOUT = downcall("wgpuBindGroupLayoutRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_PIPELINE_LAYOUT = downcall("wgpuPipelineLayoutRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_COMPUTE_PIPELINE = downcall("wgpuComputePipelineRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_SHADER_MODULE = downcall("wgpuShaderModuleRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_BIND_GROUP = downcall("wgpuBindGroupRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_COMMAND_ENCODER = downcall("wgpuCommandEncoderRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_COMPUTE_PASS_ENCODER = downcall("wgpuComputePassEncoderRelease",
        FunctionDescriptor.ofVoid(PTR));
    private static final MethodHandle MH_RELEASE_COMMAND_BUFFER = downcall("wgpuCommandBufferRelease",
        FunctionDescriptor.ofVoid(PTR));

    private static final MemorySegment ADAPTER_CALLBACK_STUB = upcall("adapterCallback",
        FunctionDescriptor.ofVoid(I32, PTR, STRING_VIEW, PTR, PTR));
    private static final MemorySegment DEVICE_CALLBACK_STUB = upcall("deviceCallback",
        FunctionDescriptor.ofVoid(I32, PTR, STRING_VIEW, PTR, PTR));
    private static final MemorySegment MAP_CALLBACK_STUB = upcall("mapCallback",
        FunctionDescriptor.ofVoid(I32, STRING_VIEW, PTR, PTR));

    private static final AtomicLong TOKEN_SEQ = new AtomicLong();
    private static final ConcurrentHashMap<Long, Pending> PENDING = new ConcurrentHashMap<>();

    private static volatile DeviceContext CTX;

    private WebGpuPanama() {}

    static double[] computeObjectiveVals(double[] WFlat, double[] KFlat, int r, int n) {
        DeviceContext ctx = context();
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment w = toNativeFloatSegment(arena, WFlat);
            MemorySegment k = toNativeFloatSegment(arena, KFlat);
            MemorySegment result = runObjectivePipeline(ctx, arena, w, k, r, n);
            return readFloatSegment(result, r);
        } catch (Throwable t) {
            throw rethrow(t);
        }
    }

    static double[] computeRandomizationMetrics(int[] designsFlat, int r, int n) {
        DeviceContext ctx = context();
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment designs = toNativeIntSegment(arena, designsFlat);
            MemorySegment result = runRandomizationPipeline(ctx, arena, designs, r, n);
            return readFloatSegment(result, n * n);
        } catch (Throwable t) {
            throw rethrow(t);
        }
    }

    public static class GpuGreedySearchSession implements AutoCloseable {
        private final DeviceContext ctx;
        private final Arena arena;
        private final MemorySegment xBuf, sinvBuf, iTsBuf, iCsBuf, avgTBuf, avgCBuf, resultsBuf, stagingBuf, paramsBuf;
        private final MemorySegment pipeline, bindGroup;
        private final int n, p, numIT, numIC;

        GpuGreedySearchSession(DeviceContext ctx, int n, int p, double[] XFlat, double[] SinvFlat, int[] i_Ts, int[] i_Cs) throws Throwable {
            this.ctx = ctx;
            this.arena = Arena.ofConfined();
            this.n = n; this.p = p;
            this.numIT = i_Ts.length; this.numIC = i_Cs.length;

            long xSize = XFlat.length * 4L, sinvSize = SinvFlat.length * 4L;
            long iSize = (long)numIT * 4, cSize = (long)numIC * 4, pSize = (long)p * 4;
            long resSize = (long)numIT * numIC * 4;

            xBuf = createBuffer(ctx.device, arena, xSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "gs_x");
            sinvBuf = createBuffer(ctx.device, arena, sinvSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "gs_sinv");
            iTsBuf = createBuffer(ctx.device, arena, iSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "gs_iTs");
            iCsBuf = createBuffer(ctx.device, arena, cSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "gs_iCs");
            avgTBuf = createBuffer(ctx.device, arena, pSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "gs_avgT");
            avgCBuf = createBuffer(ctx.device, arena, pSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "gs_avgC");
            resultsBuf = createBuffer(ctx.device, arena, resSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, "gs_res");
            stagingBuf = createBuffer(ctx.device, arena, resSize, WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst, "gs_staging");
            paramsBuf = createBuffer(ctx.device, arena, 32L, WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst, "gs_params");

            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, xBuf, 0L, toNativeFloatSegment(arena, XFlat), xSize);
            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, sinvBuf, 0L, toNativeFloatSegment(arena, SinvFlat), sinvSize);
            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, iTsBuf, 0L, toNativeIntSegment(arena, i_Ts), iSize);
            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, iCsBuf, 0L, toNativeIntSegment(arena, i_Cs), cSize);

            MemorySegment shader = createShaderModule(ctx.device, arena, "greedy", GpuBridgeResource.greedySearch());
            MemorySegment bgl = createGreedyBindGroupLayout(ctx.device, arena);
            MemorySegment pl = createPipelineLayout(ctx.device, arena, bgl);
            pipeline = createComputePipeline(ctx.device, arena, pl, shader, "main");
            bindGroup = createGreedyBindGroup(ctx.device, arena, bgl, paramsBuf, xBuf, sinvBuf, iTsBuf, iCsBuf, avgTBuf, avgCBuf, resultsBuf);
        }

        public double[] runIteration(int nT, double[] currentAvgT, double[] currentAvgC) throws Throwable {
            MemoryLayout paramsLayout = MemoryLayout.structLayout(I32.withName("n"), I32.withName("p"), I32.withName("nT"), I32.withName("num_i_Ts"), I32.withName("num_i_Cs"), MemoryLayout.paddingLayout(4), MemoryLayout.paddingLayout(8));
            MemorySegment pSeg = zero(arena, paramsLayout);
            setI32(pSeg, paramsLayout, "n", n); setI32(pSeg, paramsLayout, "p", p); setI32(pSeg, paramsLayout, "nT", nT);
            setI32(pSeg, paramsLayout, "num_i_Ts", numIT); setI32(pSeg, paramsLayout, "num_i_Cs", numIC);
            
            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, paramsBuf, 0L, pSeg, pSeg.byteSize());
            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, avgTBuf, 0L, toNativeFloatSegment(arena, currentAvgT), (long)p * 4);
            MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, avgCBuf, 0L, toNativeFloatSegment(arena, currentAvgC), (long)p * 4);

            long resSize = (long)numIT * numIC * 4;
            long workgroups = (resSize / 4L + 63L) / 64L;
            MemorySegment submission = dispatchAndCopy(ctx, arena, pipeline, bindGroup, resultsBuf, stagingBuf, workgroups, resSize);
            MemorySegment result = mapAndRead(ctx, arena, stagingBuf, numIT * numIC, submission);
            
            double[] out = new double[numIT * numIC];
            for (int i = 0; i < out.length; i++) out[i] = result.get(ValueLayout.JAVA_FLOAT, i * 4L);
            MH_BUFFER_UNMAP.invoke(stagingBuf);
            return out;
        }

        @Override
        public void close() {
            arena.close();
        }
    }

    public static GpuGreedySearchSession createGreedySearchSession(int n, int p, double[] XFlat, double[] SinvFlat, int[] i_Ts, int[] i_Cs) {
        try {
            return new GpuGreedySearchSession(context(), n, p, XFlat, SinvFlat, i_Ts, i_Cs);
        } catch (Throwable t) {
            throw rethrow(t);
        }
    }

    private static MemorySegment createGreedyBindGroupLayout(MemorySegment device, Arena arena) throws Throwable {
        MemoryLayout entryLayout = bindGroupLayoutEntry();
        MemoryLayout layoutLayout = bindGroupLayoutDescriptor();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 8, entryLayout.byteAlignment());
        for (int i = 0; i < 8; i++) {
            MemorySegment e = entries.asSlice((long) i * entryLayout.byteSize(), entryLayout.byteSize());
            setI32(e, entryLayout, "binding", i);
            setLong(e, entryLayout, "visibility", WGPUShaderStage_Compute);
            setI32(e, entryLayout, "buffer_type", i == 0 ? WGPUBufferBindingType_Uniform : (i == 7 ? WGPUBufferBindingType_Storage : WGPUBufferBindingType_ReadOnlyStorage));
        }
        MemorySegment desc = zero(arena, layoutLayout);
        setStringViewNullFlat(desc, layoutLayout, "label", arena);
        setLong(desc, layoutLayout, "entryCount", 8L);
        setAddr(desc, layoutLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP_LAYOUT.invoke(device, desc);
    }

    private static MemorySegment createGreedyBindGroup(MemorySegment device, Arena arena, MemorySegment layout, MemorySegment params, MemorySegment x, MemorySegment sinv, MemorySegment iTs, MemorySegment iCs, MemorySegment avgT, MemorySegment avgC, MemorySegment res) throws Throwable {
        MemoryLayout entryLayout = bindGroupEntryLayout();
        MemoryLayout descriptorLayout = bindGroupDescriptorLayout();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 8, entryLayout.byteAlignment());
        MemorySegment[] bufs = {params, x, sinv, iTs, iCs, avgT, avgC, res};
        for (int i = 0; i < 8; i++) {
            fillBindGroupEntry(entries.asSlice((long) i * entryLayout.byteSize()), i, bufs[i], 0L, WGPU_WHOLE_SIZE, MemorySegment.NULL, MemorySegment.NULL);
        }
        MemorySegment desc = zero(arena, descriptorLayout);
        setAddr(desc, descriptorLayout, "layout", layout);
        setLong(desc, descriptorLayout, "entryCount", 8L);
        setAddr(desc, descriptorLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP.invoke(device, desc);
    }

    static double[] computeGreedySwapObjectives(double[] XFlat, double[] SinvFlat, int n, int p, int nT, int[] indicT) {
        DeviceContext ctx = context();
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment x = toNativeFloatSegment(arena, XFlat);
            MemorySegment sinv = toNativeFloatSegment(arena, SinvFlat);
            MemorySegment result = runGreedyPipeline(ctx, arena, x, sinv, n, p, nT, toNativeIntSegment(arena, indicT));
            int numSwaps = countSwaps(indicT);
            return readFloatSegment(result, numSwaps);
        } catch (Throwable t) {
            throw rethrow(t);
        }
    }

    private static MemorySegment runGreedyPipeline(DeviceContext ctx, Arena arena, MemorySegment x, MemorySegment sinv, int n, int p, int nT, MemorySegment indicT) throws Throwable {
        // This is a stub for the full greedy pipeline if needed outside the session
        throw new UnsupportedOperationException("Use GpuGreedySearchSession for greedy search");
    }

    static double[] computeOptimalObjectives(double[] XFlat, double[] SinvFlat, int n, int p) {
        DeviceContext ctx = context();
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment x = toNativeFloatSegment(arena, XFlat);
            MemorySegment sinv = toNativeFloatSegment(arena, SinvFlat);
            MemorySegment result = runOptimalPipeline(ctx, arena, x, sinv, n, p);
            long numDesigns = choose(n, n / 2);
            if (numDesigns > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("Too many optimal designs to materialize: " + numDesigns);
            }
            return readFloatSegment(result, (int) numDesigns);
        } catch (Throwable t) {
            throw rethrow(t);
        }
    }

    private static MemorySegment runOptimalPipeline(DeviceContext ctx, Arena arena, MemorySegment x, MemorySegment sinv, int n, int p) throws Throwable {
        throw new UnsupportedOperationException("Optimal search pipeline not implemented in Java-GPU-Panama yet");
    }

    public static double[] computeMultipleKernelSwapObjectives(int n, int m, double maximum_gain_scaling, double[] packedKgrams, double[] weights, double[] initialObjs, double[] runningSums, double[] maxReds, int[] w, int[] i_Ts, int[] i_Cs) {
        DeviceContext ctx = context();
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment kgramsSeg = toNativeFloatSegment(arena, packedKgrams);
            MemorySegment weightsSeg = toNativeFloatSegment(arena, weights);
            MemorySegment initialObjsSeg = toNativeFloatSegment(arena, initialObjs);
            MemorySegment runningSumsSeg = toNativeFloatSegment(arena, runningSums);
            MemorySegment maxRedsSeg = toNativeFloatSegment(arena, maxReds);
            MemorySegment wSeg = toNativeFloatSegment(arena, toDoubleArray(w));
            MemorySegment iTsSeg = toNativeIntSegment(arena, i_Ts);
            MemorySegment iCsSeg = toNativeIntSegment(arena, i_Cs);
            
            MemorySegment result = runMultipleKernelPipeline(ctx, arena, n, m, maximum_gain_scaling, kgramsSeg, weightsSeg, initialObjsSeg, runningSumsSeg, maxRedsSeg, wSeg, iTsSeg, iCsSeg);
            return readFloatSegment(result, i_Ts.length * i_Cs.length);
        } catch (Throwable t) {
            throw rethrow(t);
        }
    }

    private static MemorySegment runMultipleKernelPipeline(DeviceContext ctx, Arena arena, int n, int m, double scaling, MemorySegment kgrams, MemorySegment weights, MemorySegment initialObjs, MemorySegment runningSums, MemorySegment maxReds, MemorySegment w, MemorySegment iTs, MemorySegment iCs) throws Throwable {
        MemorySegment shader = createShaderModule(ctx.device, arena, "multiple_kernel", GpuBridgeResource.multipleKernelSearch());
        MemoryLayout paramsLayout = MemoryLayout.structLayout(I32.withName("n"), I32.withName("m"), I32.withName("num_i_Ts"), I32.withName("num_i_Cs"), ValueLayout.JAVA_FLOAT.withName("maximum_gain_scaling"), MemoryLayout.paddingLayout(4), MemoryLayout.paddingLayout(8));
        MemorySegment params = zero(arena, paramsLayout);
        setI32(params, paramsLayout, "n", n);
        setI32(params, paramsLayout, "m", m);
        setI32(params, paramsLayout, "num_i_Ts", (int) iTs.byteSize() / 4);
        setI32(params, paramsLayout, "num_i_Cs", (int) iCs.byteSize() / 4);
        params.set(ValueLayout.JAVA_FLOAT, off(paramsLayout, "maximum_gain_scaling"), (float) scaling);
        
        long resSize = (long) ((int) iTs.byteSize() / 4) * ((int) iCs.byteSize() / 4) * 4L;
        MemorySegment paramsBuf = createBuffer(ctx.device, arena, params.byteSize(), WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst, "mk_params");
        MemorySegment kBuf = createBuffer(ctx.device, arena, kgrams.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_kgrams");
        MemorySegment weightsBuf = createBuffer(ctx.device, arena, weights.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_weights");
        MemorySegment initBuf = createBuffer(ctx.device, arena, initialObjs.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_init");
        MemorySegment runBuf = createBuffer(ctx.device, arena, runningSums.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_run");
        MemorySegment maxRedBuf = createBuffer(ctx.device, arena, maxReds.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_maxred");
        MemorySegment wBuf = createBuffer(ctx.device, arena, w.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_w");
        MemorySegment iTsBuf = createBuffer(ctx.device, arena, iTs.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_iTs");
        MemorySegment iCsBuf = createBuffer(ctx.device, arena, iCs.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "mk_iCs");
        MemorySegment outBuf = createBuffer(ctx.device, arena, resSize, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, "mk_out");
        MemorySegment staging = createBuffer(ctx.device, arena, resSize, WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst, "mk_staging");

        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, paramsBuf, 0L, params, params.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, kBuf, 0L, kgrams, kgrams.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, weightsBuf, 0L, weights, weights.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, initBuf, 0L, initialObjs, initialObjs.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, runBuf, 0L, runningSums, runningSums.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, maxRedBuf, 0L, maxReds, maxReds.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, wBuf, 0L, w, w.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, iTsBuf, 0L, iTs, iTs.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, iCsBuf, 0L, iCs, iCs.byteSize());

        MemorySegment bgl = createMultipleKernelBindGroupLayout(ctx.device, arena);
        MemorySegment pl = createPipelineLayout(ctx.device, arena, bgl);
        MemorySegment pipeline = createComputePipeline(ctx.device, arena, pl, shader, "main");
        MemorySegment bindGroup = createMultipleKernelBindGroup(ctx.device, arena, bgl, paramsBuf, kBuf, weightsBuf, initBuf, runBuf, maxRedBuf, wBuf, iTsBuf, iCsBuf, outBuf);
        
        long workgroups = (resSize / 4L + 63L) / 64L;
        MemorySegment submission = dispatchAndCopy(ctx, arena, pipeline, bindGroup, outBuf, staging, workgroups, resSize);
        MemorySegment result = mapAndRead(ctx, arena, staging, (int) (resSize / 4L), submission);
        
        return result;
    }

    private static MemorySegment createMultipleKernelBindGroupLayout(MemorySegment device, Arena arena) throws Throwable {
        MemoryLayout entryLayout = bindGroupLayoutEntry();
        MemoryLayout layoutLayout = bindGroupLayoutDescriptor();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 10, entryLayout.byteAlignment());
        for (int i = 0; i < 10; i++) {
            MemorySegment e = entries.asSlice((long) i * entryLayout.byteSize(), entryLayout.byteSize());
            setI32(e, entryLayout, "binding", i);
            setLong(e, entryLayout, "visibility", WGPUShaderStage_Compute);
            setI32(e, entryLayout, "buffer_type", i == 0 ? WGPUBufferBindingType_Uniform : (i == 9 ? WGPUBufferBindingType_Storage : WGPUBufferBindingType_ReadOnlyStorage));
        }
        MemorySegment desc = zero(arena, layoutLayout);
        setLong(desc, layoutLayout, "entryCount", 10L);
        setAddr(desc, layoutLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP_LAYOUT.invoke(device, desc);
    }

    private static MemorySegment createMultipleKernelBindGroup(MemorySegment device, Arena arena, MemorySegment layout, MemorySegment params, MemorySegment k, MemorySegment wts, MemorySegment init, MemorySegment run, MemorySegment maxr, MemorySegment w, MemorySegment iTs, MemorySegment iCs, MemorySegment out) throws Throwable {
        MemoryLayout entryLayout = bindGroupEntryLayout();
        MemoryLayout descriptorLayout = bindGroupDescriptorLayout();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 10, entryLayout.byteAlignment());
        MemorySegment[] bufs = {params, k, wts, init, run, maxr, w, iTs, iCs, out};
        for (int i = 0; i < 10; i++) {
            fillBindGroupEntry(entries.asSlice((long) i * entryLayout.byteSize()), i, bufs[i], 0L, WGPU_WHOLE_SIZE, MemorySegment.NULL, MemorySegment.NULL);
        }
        MemorySegment desc = zero(arena, descriptorLayout);
        setAddr(desc, descriptorLayout, "layout", layout);
        setLong(desc, descriptorLayout, "entryCount", 10L);
        setAddr(desc, descriptorLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP.invoke(device, desc);
    }

    private static double[] toDoubleArray(int[] vals) {
        double[] out = new double[vals.length];
        for (int i = 0; i < vals.length; i++) out[i] = (double) vals[i];
        return out;
    }

    public static double computeObjective(double[] WFlat, double[] KFlat, int r, int n) {
        double[] vals = computeObjectiveVals(WFlat, KFlat, r, n);
        return vals.length == 0 ? Double.NaN : vals[0];
    }

    private static void adapterCallback(int status, MemorySegment adapter, MemorySegment message, MemorySegment userdata1, MemorySegment userdata2) {
        completePending(status, adapter, userdata1);
    }

    private static void deviceCallback(int status, MemorySegment device, MemorySegment message, MemorySegment userdata1, MemorySegment userdata2) {
        completePending(status, device, userdata1);
    }

    private static void mapCallback(int status, MemorySegment message, MemorySegment userdata1, MemorySegment userdata2) {
        completePending(status, MemorySegment.NULL, userdata1);
    }

    private static void completePending(int status, MemorySegment handle, MemorySegment userdata1) {
        long token = MemorySegment.ofAddress(userdata1.address()).reinterpret(I64.byteSize()).get(I64, 0);
        Pending pending = PENDING.remove(token);
        if (pending != null) {
            pending.status = status;
            pending.handle = handle;
            pending.latch.countDown();
        }
    }

    private static synchronized DeviceContext context() {
        if (CTX == null) {
            CTX = initContext();
        }
        return CTX;
    }

    private static DeviceContext initContext() {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment instanceDesc = zero(arena, WGPU_INSTANCE_DESCRIPTOR);
            MemorySegment instance = (MemorySegment) MH_CREATE_INSTANCE.invoke(instanceDesc);
            if (instance.equals(MemorySegment.NULL)) {
                throw new IllegalStateException("wgpuCreateInstance returned null");
            }

            Pending adapterPending = new Pending();
            MemorySegment adapterToken = registerPending(arena, adapterPending);
            MemorySegment adapterOpts = zero(arena, REQUEST_ADAPTER_OPTIONS);
            setI32(adapterOpts, REQUEST_ADAPTER_OPTIONS, "featureLevel", WGPUFeatureLevel_Undefined);
            setI32(adapterOpts, REQUEST_ADAPTER_OPTIONS, "powerPreference", WGPUPowerPreference_HighPerformance);
            setI32(adapterOpts, REQUEST_ADAPTER_OPTIONS, "forceFallbackAdapter", WGPU_FALSE);
            setI32(adapterOpts, REQUEST_ADAPTER_OPTIONS, "backendType", WGPUBackendType_Undefined);
            setAddr(adapterOpts, REQUEST_ADAPTER_OPTIONS, "compatibleSurface", MemorySegment.NULL);
            MemorySegment adapterCb = zero(arena, REQUEST_ADAPTER_CALLBACK_INFO);
            setI32(adapterCb, REQUEST_ADAPTER_CALLBACK_INFO, "mode", WGPUCallbackMode_AllowSpontaneous);
            setAddr(adapterCb, REQUEST_ADAPTER_CALLBACK_INFO, "callback", ADAPTER_CALLBACK_STUB);
            setAddr(adapterCb, REQUEST_ADAPTER_CALLBACK_INFO, "userdata1", adapterToken);
            setAddr(adapterCb, REQUEST_ADAPTER_CALLBACK_INFO, "userdata2", MemorySegment.NULL);
            MH_INSTANCE_REQUEST_ADAPTER.invoke(instance, adapterOpts, adapterCb);
            await(adapterPending);
            if (adapterPending.status != WGPURequestAdapterStatus_Success || adapterPending.handle.equals(MemorySegment.NULL)) {
                throw new IllegalStateException("No usable WebGPU adapter available");
            }

            Pending devicePending = new Pending();
            MemorySegment deviceToken = registerPending(arena, devicePending);
            MemorySegment deviceDesc = zero(arena, DEVICE_DESCRIPTOR);
            setStringViewFlat(deviceDesc, DEVICE_DESCRIPTOR, "label", arena, "GreedyExperimentalDesign WebGPU device");
            setAddr(deviceDesc, DEVICE_DESCRIPTOR, "defaultQueue_nextInChain", MemorySegment.NULL);
            setStringViewNullFlat(deviceDesc, DEVICE_DESCRIPTOR, "defaultQueue_label", arena);
            MemorySegment deviceCb = zero(arena, REQUEST_DEVICE_CALLBACK_INFO);
            setI32(deviceCb, REQUEST_DEVICE_CALLBACK_INFO, "mode", WGPUCallbackMode_AllowSpontaneous);
            setAddr(deviceCb, REQUEST_DEVICE_CALLBACK_INFO, "callback", DEVICE_CALLBACK_STUB);
            setAddr(deviceCb, REQUEST_DEVICE_CALLBACK_INFO, "userdata1", deviceToken);
            setAddr(deviceCb, REQUEST_DEVICE_CALLBACK_INFO, "userdata2", MemorySegment.NULL);
            MH_ADAPTER_REQUEST_DEVICE.invoke(adapterPending.handle, deviceDesc, deviceCb);
            await(devicePending);
            if (devicePending.status != WGPURequestDeviceStatus_Success || devicePending.handle.equals(MemorySegment.NULL)) {
                throw new IllegalStateException("No usable WebGPU device available");
            }

            MemorySegment queue = (MemorySegment) MH_DEVICE_GET_QUEUE.invoke(devicePending.handle);
            if (queue.equals(MemorySegment.NULL)) {
                throw new IllegalStateException("Failed to get WebGPU queue");
            }
            return new DeviceContext(instance, adapterPending.handle, devicePending.handle, queue);
        } catch (Throwable t) {
            throw new RuntimeException(t);
        }
    }

    private static MemorySegment runObjectivePipeline(DeviceContext ctx, Arena arena, MemorySegment w, MemorySegment k, int r, int n) throws Throwable {
        MemorySegment shader = createShaderModule(ctx.device, arena, "objective", GpuBridgeResource.kernelQuadraticForm());
        MemoryLayout paramsLayout = MemoryLayout.structLayout(I32.withName("n"), I32.withName("r"), I32.withName("pad0"), I32.withName("pad1"));
        MemorySegment params = zero(arena, paramsLayout);
        setI32(params, paramsLayout, "n", n);
        setI32(params, paramsLayout, "r", r);
        MemorySegment paramsBuf = createBuffer(ctx.device, arena, params.byteSize(), WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst, "objective_params");
        MemorySegment wBuf = createBuffer(ctx.device, arena, w.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "objective_w");
        MemorySegment kBuf = createBuffer(ctx.device, arena, k.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "objective_k");
        MemorySegment outBuf = createBuffer(ctx.device, arena, (long) r * 4L, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, "objective_out");
        MemorySegment staging = createBuffer(ctx.device, arena, (long) r * 4L, WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst, "objective_staging");
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, paramsBuf, 0L, params, params.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, wBuf, 0L, w, w.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, kBuf, 0L, k, k.byteSize());
        MemorySegment bgl = createObjectiveBindGroupLayout(ctx.device, arena, r, n);
        MemorySegment pl = createPipelineLayout(ctx.device, arena, bgl);
        MemorySegment pipeline = createComputePipeline(ctx.device, arena, pl, shader, "main");
        MemorySegment bindGroup = createObjectiveBindGroup(ctx.device, arena, bgl, paramsBuf, kBuf, wBuf, outBuf);
        MemorySegment submission = dispatchAndCopy(ctx, arena, pipeline, bindGroup, outBuf, staging, (r + 63L) / 64L, r * 4L);
        MemorySegment result = mapAndRead(ctx, arena, staging, r, submission);
        return result;
    }

    private static MemorySegment runRandomizationPipeline(DeviceContext ctx, Arena arena, MemorySegment designs, int r, int n) throws Throwable {
        MemorySegment shader = createShaderModule(ctx.device, arena, "randomization", GpuBridgeResource.randomizationMetrics());
        MemoryLayout paramsLayout = MemoryLayout.structLayout(I32.withName("n"), I32.withName("r"), I32.withName("pad0"), I32.withName("pad1"));
        MemorySegment params = zero(arena, paramsLayout);
        setI32(params, paramsLayout, "n", n);
        setI32(params, paramsLayout, "r", r);
        MemorySegment paramsBuf = createBuffer(ctx.device, arena, params.byteSize(), WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst, "randomization_params");
        MemorySegment designsBuf = createBuffer(ctx.device, arena, designs.byteSize(), WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, "randomization_designs");
        MemorySegment outBuf = createBuffer(ctx.device, arena, (long) n * (long) n * 4L, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, "randomization_out");
        MemorySegment staging = createBuffer(ctx.device, arena, (long) n * (long) n * 4L, WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst, "randomization_staging");
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, paramsBuf, 0L, params, params.byteSize());
        MH_QUEUE_WRITE_BUFFER.invoke(ctx.queue, designsBuf, 0L, designs, designs.byteSize());
        MemorySegment bgl = createRandomizationBindGroupLayout(ctx.device, arena, r, n);
        MemorySegment pl = createPipelineLayout(ctx.device, arena, bgl);
        MemorySegment pipeline = createComputePipeline(ctx.device, arena, pl, shader, "main");
        MemorySegment bindGroup = createRandomizationBindGroup(ctx.device, arena, bgl, paramsBuf, designsBuf, outBuf);
        MemorySegment submission = dispatchAndCopy(ctx, arena, pipeline, bindGroup, outBuf, staging, ((long) n * (n - 1L) / 2L + 63L) / 64L, (long) n * (long) n * 4L);
        MemorySegment result = mapAndRead(ctx, arena, staging, n * n, submission);
        return result;
    }

    private static MemorySegment createObjectiveBindGroupLayout(MemorySegment device, Arena arena, int r, int n) throws Throwable {
        MemoryLayout entryLayout = bindGroupLayoutEntry();
        MemoryLayout layoutLayout = bindGroupLayoutDescriptor();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 4, entryLayout.byteAlignment());
        fillBindGroupLayoutEntry(entries.asSlice(0), entryLayout, 0, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_Uniform, 8L);
        fillBindGroupLayoutEntry(entries.asSlice(entryLayout.byteSize()), entryLayout, 1, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_ReadOnlyStorage, (long) n * (long) n * 4L);
        fillBindGroupLayoutEntry(entries.asSlice(entryLayout.byteSize() * 2L), entryLayout, 2, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_ReadOnlyStorage, (long) r * (long) n * 4L);
        fillBindGroupLayoutEntry(entries.asSlice(entryLayout.byteSize() * 3L), entryLayout, 3, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_Storage, (long) r * 4L);
        MemorySegment desc = zero(arena, layoutLayout);
        setLong(desc, layoutLayout, "entryCount", 4L);
        setAddr(desc, layoutLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP_LAYOUT.invoke(device, desc);
    }

    private static MemorySegment createObjectiveBindGroup(MemorySegment device, Arena arena, MemorySegment layout, MemorySegment params, MemorySegment k, MemorySegment w, MemorySegment out) throws Throwable {
        MemoryLayout entryLayout = bindGroupEntryLayout();
        MemoryLayout descriptorLayout = bindGroupDescriptorLayout();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 4, entryLayout.byteAlignment());
        MemorySegment[] bufs = {params, k, w, out};
        for (int i = 0; i < 4; i++) {
            fillBindGroupEntry(entries.asSlice((long) i * entryLayout.byteSize()), i, bufs[i], 0L, WGPU_WHOLE_SIZE, MemorySegment.NULL, MemorySegment.NULL);
        }
        MemorySegment desc = zero(arena, descriptorLayout);
        setAddr(desc, descriptorLayout, "layout", layout);
        setLong(desc, descriptorLayout, "entryCount", 4L);
        setAddr(desc, descriptorLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP.invoke(device, desc);
    }

    private static MemorySegment createRandomizationBindGroupLayout(MemorySegment device, Arena arena, int r, int n) throws Throwable {
        MemoryLayout entryLayout = bindGroupLayoutEntry();
        MemoryLayout layoutLayout = bindGroupLayoutDescriptor();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 3, entryLayout.byteAlignment());
        fillBindGroupLayoutEntry(entries.asSlice(0), entryLayout, 0, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_Uniform, 8L);
        fillBindGroupLayoutEntry(entries.asSlice(entryLayout.byteSize()), entryLayout, 1, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_ReadOnlyStorage, (long) r * (long) n * 4L);
        fillBindGroupLayoutEntry(entries.asSlice(entryLayout.byteSize() * 2L), entryLayout, 2, (int) WGPUShaderStage_Compute, WGPUBufferBindingType_Storage, (long) n * (long) n * 4L);
        MemorySegment desc = zero(arena, layoutLayout);
        setLong(desc, layoutLayout, "entryCount", 3L);
        setAddr(desc, layoutLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP_LAYOUT.invoke(device, desc);
    }

    private static MemorySegment createRandomizationBindGroup(MemorySegment device, Arena arena, MemorySegment layout, MemorySegment params, MemorySegment designs, MemorySegment out) throws Throwable {
        MemoryLayout entryLayout = bindGroupEntryLayout();
        MemoryLayout descriptorLayout = bindGroupDescriptorLayout();
        MemorySegment entries = zero(arena, entryLayout.byteSize() * 3, entryLayout.byteAlignment());
        MemorySegment[] bufs = {params, designs, out};
        for (int i = 0; i < 3; i++) {
            fillBindGroupEntry(entries.asSlice((long) i * entryLayout.byteSize()), i, bufs[i], 0L, WGPU_WHOLE_SIZE, MemorySegment.NULL, MemorySegment.NULL);
        }
        MemorySegment desc = zero(arena, descriptorLayout);
        setAddr(desc, descriptorLayout, "layout", layout);
        setLong(desc, descriptorLayout, "entryCount", 3L);
        setAddr(desc, descriptorLayout, "entries", entries);
        return (MemorySegment) MH_DEVICE_CREATE_BIND_GROUP.invoke(device, desc);
    }

    private static MemorySegment createShaderModule(MemorySegment device, Arena arena, String label, String code) throws Throwable {
        MemorySegment codeSeg = arena.allocateFrom(code);
        MemorySegment wgsl = zero(arena, SHADER_SOURCE_WGSL);
        setI32(wgsl, SHADER_SOURCE_WGSL, "sType", 2); // WGPUSType_ShaderSourceWGSL = 0x00000002
        setAddr(wgsl, SHADER_SOURCE_WGSL, "code_data", codeSeg);
        setLong(wgsl, SHADER_SOURCE_WGSL, "code_length", WGPU_STRLEN);
        MemorySegment desc = zero(arena, SHADER_MODULE_DESCRIPTOR);
        setAddr(desc, SHADER_MODULE_DESCRIPTOR, "nextInChain", wgsl);
        setStringViewFlat(desc, SHADER_MODULE_DESCRIPTOR, "label", arena, label);
        return (MemorySegment) MH_DEVICE_CREATE_SHADER_MODULE.invoke(device, desc);
    }

    private static MemorySegment createBuffer(MemorySegment device, Arena arena, long size, long usage, String label) throws Throwable {
        MemorySegment desc = zero(arena, BUFFER_DESCRIPTOR);
        setStringViewFlat(desc, BUFFER_DESCRIPTOR, "label", arena, label);
        setLong(desc, BUFFER_DESCRIPTOR, "usage", usage);
        setLong(desc, BUFFER_DESCRIPTOR, "size", size);
        setI32(desc, BUFFER_DESCRIPTOR, "mappedAtCreation", WGPU_FALSE);
        return (MemorySegment) MH_DEVICE_CREATE_BUFFER.invoke(device, desc);
    }

    private static MemoryLayout bindGroupLayoutEntry() {
        // Matches WGPUBindGroupLayoutEntry from webgpu.h (120 bytes on 64-bit).
        // visibility is WGPUShaderStage = WGPUFlags = uint64_t.
        // bindingArraySize is a field between visibility and the buffer sub-struct.
        return MemoryLayout.structLayout(
            PTR.withName("nextInChain"),
            I32.withName("binding"),
            MemoryLayout.paddingLayout(4),
            I64.withName("visibility"),
            I32.withName("bindingArraySize"),
            MemoryLayout.paddingLayout(4),
            PTR.withName("buffer_nextInChain"),
            I32.withName("buffer_type"),
            I32.withName("buffer_hasDynamicOffset"),
            I64.withName("buffer_minBindingSize"),
            PTR.withName("sampler_nextInChain"),
            I32.withName("sampler_type"),
            MemoryLayout.paddingLayout(4),
            PTR.withName("texture_nextInChain"),
            I32.withName("texture_sampleType"),
            I32.withName("texture_viewDimension"),
            I32.withName("texture_multisampled"),
            MemoryLayout.paddingLayout(4),
            PTR.withName("storageTexture_nextInChain"),
            I32.withName("storageTexture_access"),
            I32.withName("storageTexture_format"),
            I32.withName("storageTexture_viewDimension"),
            MemoryLayout.paddingLayout(4)
        );
    }

    private static MemoryLayout bindGroupLayoutDescriptor() {
        return MemoryLayout.structLayout(PTR.withName("nextInChain"), PTR.withName("label_data"), I64.withName("label_length"), I64.withName("entryCount"), PTR.withName("entries"));
    }

    private static MemoryLayout bindGroupEntryLayout() {
        return MemoryLayout.structLayout(PTR.withName("nextInChain"), I32.withName("binding"), MemoryLayout.paddingLayout(4), PTR.withName("buffer"), I64.withName("offset"), I64.withName("size"), PTR.withName("sampler"), PTR.withName("textureView"));
    }

    private static MemoryLayout bindGroupDescriptorLayout() {
        return MemoryLayout.structLayout(PTR.withName("nextInChain"), PTR.withName("label_data"), I64.withName("label_length"), PTR.withName("layout"), I64.withName("entryCount"), PTR.withName("entries"));
    }

    private static MemorySegment createPipelineLayout(MemorySegment device, Arena arena, MemorySegment bindGroupLayout) throws Throwable {
        MemorySegment layouts = arena.allocate(PTR.byteSize(), PTR.byteAlignment()).fill((byte) 0);
        layouts.set(PTR, 0, bindGroupLayout);
        MemorySegment desc = zero(arena, PIPELINE_LAYOUT_DESCRIPTOR);
        setLong(desc, PIPELINE_LAYOUT_DESCRIPTOR, "bindGroupLayoutCount", 1L);
        setAddr(desc, PIPELINE_LAYOUT_DESCRIPTOR, "bindGroupLayouts", layouts);
        return (MemorySegment) MH_DEVICE_CREATE_PIPELINE_LAYOUT.invoke(device, desc);
    }

    private static MemorySegment createComputePipeline(MemorySegment device, Arena arena, MemorySegment pipelineLayout, MemorySegment shader, String entryPoint) throws Throwable {
        MemorySegment desc = zero(arena, COMPUTE_PIPELINE_DESCRIPTOR);
        setAddr(desc, COMPUTE_PIPELINE_DESCRIPTOR, "layout", pipelineLayout);
        setAddr(desc, COMPUTE_PIPELINE_DESCRIPTOR, "compute_module", shader);
        setStringViewFlat(desc, COMPUTE_PIPELINE_DESCRIPTOR, "compute_entryPoint", arena, entryPoint);
        return (MemorySegment) MH_DEVICE_CREATE_COMPUTE_PIPELINE.invoke(device, desc);
    }

    private static MemorySegment dispatchAndCopy(DeviceContext ctx, Arena arena, MemorySegment pipeline, MemorySegment bindGroup, MemorySegment outBuf, MemorySegment staging, long workgroups, long copyBytes) throws Throwable {
        MemorySegment encoder = (MemorySegment) MH_DEVICE_CREATE_COMMAND_ENCODER.invoke(ctx.device, zero(arena, COMMAND_ENCODER_DESCRIPTOR));
        MemorySegment pass = (MemorySegment) MH_COMMAND_ENCODER_BEGIN_COMPUTE_PASS.invoke(encoder, zero(arena, COMPUTE_PASS_DESCRIPTOR));
        MH_COMPUTE_PASS_SET_PIPELINE.invoke(pass, pipeline);
        MH_COMPUTE_PASS_SET_BIND_GROUP.invoke(pass, 0, bindGroup, 0L, MemorySegment.NULL);
        MH_COMPUTE_PASS_DISPATCH.invoke(pass, (int) workgroups, 1, 1);
        MH_COMPUTE_PASS_END.invoke(pass);
        MH_COMMAND_COPY_BUFFER_TO_BUFFER.invoke(encoder, outBuf, 0L, staging, 0L, copyBytes);
        MemorySegment cmd = (MemorySegment) MH_COMMAND_ENCODER_FINISH.invoke(encoder, zero(arena, COMMAND_BUFFER_DESCRIPTOR));
        MemorySegment cmdArray = submitArray(arena, cmd);
        long submissionIndex = (long) MH_QUEUE_SUBMIT_FOR_INDEX.invoke(ctx.queue, 1L, cmdArray);
        return submissionIndexSegment(arena, submissionIndex);
    }

    private static MemorySegment mapAndRead(DeviceContext ctx, Arena arena, MemorySegment staging, int count, MemorySegment submission) throws Throwable {
        Pending pending = new Pending();
        MemorySegment token = registerPending(arena, pending);
        MemorySegment cb = zero(arena, BUFFER_MAP_CALLBACK_INFO);
        setI32(cb, BUFFER_MAP_CALLBACK_INFO, "mode", WGPUCallbackMode_AllowSpontaneous);
        setAddr(cb, BUFFER_MAP_CALLBACK_INFO, "callback", MAP_CALLBACK_STUB);
        setAddr(cb, BUFFER_MAP_CALLBACK_INFO, "userdata1", token);
        long mapBytes = count * 4L;
        MH_BUFFER_MAP_ASYNC.invoke(staging, WGPUMapMode_Read, 0L, mapBytes, cb);
        while (pending.latch.getCount() > 0) {
            MH_DEVICE_POLL.invoke(ctx.device, WGPU_TRUE, submission);
            Thread.sleep(1L);
        }
        if (pending.status != WGPUMapAsyncStatus_Success) {
            throw new IllegalStateException("wgpuBufferMapAsync failed with status " + pending.status);
        }
        MemorySegment mapped = (MemorySegment) MH_BUFFER_GET_CONST_MAPPED_RANGE.invoke(staging, 0L, mapBytes);
        return mapped.reinterpret(mapBytes);
    }

    private static MemorySegment submissionIndexSegment(Arena arena, long submissionIndex) {
        MemorySegment index = arena.allocate(I64);
        index.set(I64, 0, submissionIndex);
        return index;
    }

    private static MemorySegment submitArray(Arena arena, MemorySegment cmd) {
        MemorySegment arr = arena.allocate(PTR.byteSize(), PTR.byteAlignment()).fill((byte) 0);
        arr.set(PTR, 0, cmd);
        return arr;
    }

    private static void fillBindGroupLayoutEntry(MemorySegment entry, MemoryLayout layout, int binding, int visibility, int bufferType, long minBindingSize) {
        setI32(entry, layout, "binding", binding);
        setLong(entry, layout, "visibility", (long) visibility);
        setI32(entry, layout, "buffer_type", bufferType);
        setLong(entry, layout, "buffer_minBindingSize", minBindingSize);
    }

    private static void fillBindGroupEntry(MemorySegment entry, int binding, MemorySegment buffer, long offset, long size, MemorySegment sampler, MemorySegment textureView) {
        MemoryLayout layout = bindGroupEntryLayout();
        setI32(entry, layout, "binding", binding);
        setAddr(entry, layout, "buffer", buffer);
        setLong(entry, layout, "offset", offset);
        setLong(entry, layout, "size", size);
        setAddr(entry, layout, "sampler", sampler);
        setAddr(entry, layout, "textureView", textureView);
    }

    private static MemorySegment zero(Arena arena, MemoryLayout layout) {
        return arena.allocate(layout).fill((byte) 0);
    }

    private static MemorySegment zero(Arena arena, long size, long alignment) {
        return arena.allocate(size, alignment).fill((byte) 0);
    }

    private static void setI32(MemorySegment segment, MemoryLayout layout, String field, int value) {
        segment.set(I32, off(layout, field), value);
    }

    private static void setLong(MemorySegment segment, MemoryLayout layout, String field, long value) {
        segment.set(I64, off(layout, field), value);
    }

    private static void setAddr(MemorySegment segment, MemoryLayout layout, String field, MemorySegment value) {
        segment.set(PTR, off(layout, field), value);
    }

    private static long off(MemoryLayout layout, String... path) {
        PathElement[] elements = Arrays.stream(path).map(PathElement::groupElement).toArray(PathElement[]::new);
        return layout.byteOffset(elements);
    }

    private static MemorySegment registerPending(Arena arena, Pending pending) {
        long token = TOKEN_SEQ.incrementAndGet();
        PENDING.put(token, pending);
        MemorySegment tokenSeg = arena.allocate(I64);
        tokenSeg.set(I64, 0, token);
        return tokenSeg;
    }

    private static void await(Pending pending) {
        try {
            boolean ok = pending.latch.await(60, TimeUnit.SECONDS);
            if (!ok) {
                throw new RuntimeException("Timed out waiting for WebGPU callback");
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }
    }

    private static double[] readFloatSegment(MemorySegment segment, int count) {
        float[] floats = segment.reinterpret((long) count * 4L).toArray(ValueLayout.JAVA_FLOAT);
        double[] out = new double[floats.length];
        for (int i = 0; i < floats.length; i++) {
            out[i] = floats[i];
        }
        return out;
    }

    private static float[] toFloatArray(double[] values) {
        float[] out = new float[values.length];
        for (int i = 0; i < values.length; i++) {
            out[i] = (float) values[i];
        }
        return out;
    }

    private static MemorySegment toNativeFloatSegment(Arena arena, double[] values) {
        float[] floats = toFloatArray(values);
        MemorySegment out = arena.allocate((long) floats.length * Float.BYTES, Float.BYTES);
        for (int i = 0; i < floats.length; i++) {
            out.setAtIndex(ValueLayout.JAVA_FLOAT, i, floats[i]);
        }
        return out;
    }

    private static MemorySegment toNativeIntSegment(Arena arena, int[] values) {
        MemorySegment out = arena.allocate((long) values.length * Integer.BYTES, Integer.BYTES);
        for (int i = 0; i < values.length; i++) {
            out.setAtIndex(ValueLayout.JAVA_INT, i, values[i]);
        }
        return out;
    }

    private static int countSwaps(int[] indicT) {
        int countT = 0;
        int countC = 0;
        for (int v : indicT) {
            if (v == 1) {
                countT++;
            } else {
                countC++;
            }
        }
        return countT * countC;
    }

    private static long choose(int n, int kIn) {
        if (kIn > n) {
            return 0L;
        }
        int k = Math.min(kIn, n - kIn);
        long result = 1L;
        for (int i = 1; i <= k; i++) {
            result = (result * (n - k + i)) / i;
        }
        return result;
    }

    private static void setStringViewFlat(MemorySegment segment, MemoryLayout layout, String field, Arena arena, String s) {
        MemorySegment str = arena.allocateFrom(s);
        setAddr(segment, layout, field + "_data", str);
        setLong(segment, layout, field + "_length", WGPU_STRLEN);
    }

    private static void setStringViewNullFlat(MemorySegment segment, MemoryLayout layout, String field, Arena arena) {
        setAddr(segment, layout, field + "_data", MemorySegment.NULL);
        setLong(segment, layout, field + "_length", WGPU_STRLEN);
    }

    private static RuntimeException rethrow(Throwable t) {
        if (t instanceof RuntimeException re) {
            return re;
        }
        return new RuntimeException(t);
    }

    private static MethodHandle downcall(String symbol, FunctionDescriptor descriptor) {
        return LINKER.downcallHandle(LOOKUP.findOrThrow(symbol), descriptor);
    }

    private static MemorySegment upcall(String method, FunctionDescriptor descriptor) {
        try {
            MethodHandle handle = MethodHandles.lookup().findStatic(WebGpuPanama.class, method,
                switch (method) {
                    case "adapterCallback", "deviceCallback" -> MethodType.methodType(void.class, int.class, MemorySegment.class, MemorySegment.class, MemorySegment.class, MemorySegment.class);
                    case "mapCallback" -> MethodType.methodType(void.class, int.class, MemorySegment.class, MemorySegment.class, MemorySegment.class);
                    default -> throw new IllegalArgumentException(method);
                });
            return LINKER.upcallStub(handle, descriptor, UPCALL_ARENA);
        } catch (NoSuchMethodException | IllegalAccessException e) {
            throw new RuntimeException(e);
        }
    }

    private static Path findLibraryPath() {
        List<Path> candidates = List.of(
            Paths.get("lib", "wgpu", "lib", "libwgpu_native.so"),
            Paths.get("lib", "wgpu", "lib", "libwgpu_native.dylib"),
            Paths.get("lib", "wgpu", "lib", "wgpu_native.dll"),
            Paths.get("lib", "wgpu", "lib", "libwgpu_native.dll")
        );
        for (Path p : candidates) {
            if (Files.exists(p)) {
                return p.toAbsolutePath().normalize();
            }
        }
        throw new IllegalStateException("Unable to locate libwgpu_native");
    }

    private static void loadLibrary() {
        System.load(LIB_PATH.toString());
    }

    static {
        loadLibrary();
    }

    private static final class Pending {
        final CountDownLatch latch = new CountDownLatch(1);
        volatile int status;
        volatile MemorySegment handle = MemorySegment.NULL;
    }

    private static final class DeviceContext {
        final MemorySegment instance;
        final MemorySegment adapter;
        final MemorySegment device;
        final MemorySegment queue;

        DeviceContext(MemorySegment instance, MemorySegment adapter, MemorySegment device, MemorySegment queue) {
            this.instance = instance;
            this.adapter = adapter;
            this.device = device;
            this.queue = queue;
        }
    }

    private static final class GpuBridgeResource {
        static String kernelQuadraticForm() { return ResourceHolder.KERNEL; }
        static String randomizationMetrics() { return ResourceHolder.RANDOMIZATION; }
        static String greedySearch() { return ResourceHolder.GREEDY; }
        static String optimalSearch() { return ResourceHolder.OPTIMAL; }
        static String multipleKernelSearch() { return ResourceHolder.MULTIPLE_KERNEL; }

        private static final class ResourceHolder {
            static final String KERNEL = load("/ObjectiveFunctions/kernel_quadratic_form.wgsl");
            static final String RANDOMIZATION = load("/DesignMetrics/randomization_metrics.wgsl");
            static final String GREEDY = load("/GreedyExperimentalDesign/greedy_search.wgsl");
            static final String OPTIMAL = load("/OptimalExperimentalDesign/optimal_search.wgsl");
            static final String MULTIPLE_KERNEL = load("/ObjectiveFunctions/multiple_kernel_search.wgsl");

            private static String load(String path) {
                try (var in = WebGpuPanama.class.getResourceAsStream(path)) {
                    if (in == null) {
                        throw new IllegalStateException("Missing shader resource: " + path);
                    }
                    return new String(in.readAllBytes(), StandardCharsets.UTF_8);
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }
}
