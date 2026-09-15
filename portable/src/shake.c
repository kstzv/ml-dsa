#include "shake.h"

static void shake_absorb(struct shake_ctx *ctx);
static void shake_finalize(struct shake_ctx *ctx);
static void shake_squeeze(struct shake_ctx *ctx);

int shake_ctx_get_mem(struct shake_ctx *ctx, uint64_t *state, size_t size_mem)
{
	if(!ctx || !state || size_mem != 25) { return -1; }
	ctx->state = state;
	return 0;
}

void shake_ctx_init(struct shake_ctx *ctx,uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen)
{
    ctx->rate = 0;
    ctx->pos = 0;

    ctx->in = in;
    ctx->inlen = inlen;

    ctx->out = out;
    ctx->outlen = outlen;
}

void shake_ctx_zero(struct shake_ctx *ctx)
{
    if (ctx == NULL || ctx->state == NULL) { return; }

    volatile uint8_t *p = (volatile uint8_t *)ctx->state;
    size_t len = 200;
    while (len--) { *p++ = 0; }
    ctx->in = NULL;
    ctx->inlen = 0;
    ctx->out = NULL;
    ctx->outlen = 0;
    ctx->pos = 0;
    ctx->rate = 0;
}

void shake128(struct shake_ctx *ctx)
{
    ctx->rate = 168;

    if (ctx->pos != 0) 
    {
        shake_squeeze(ctx);
        return;
    }

    shake_absorb(ctx);
    shake_finalize(ctx);

    ctx->pos = 0;

    shake_squeeze(ctx);
}

void shake256(struct shake_ctx *ctx)
{
    ctx->rate = 136;

    if (ctx->pos != 0) 
    {
        shake_squeeze(ctx);
        return;
    }

    shake_absorb(ctx);
    shake_finalize(ctx);

    ctx->pos = 0;

    shake_squeeze(ctx);
}

static void shake_absorb(struct shake_ctx *ctx)
{
    uint8_t *s = (uint8_t *)ctx->state;

    while (ctx->inlen > 0)
    {
        size_t n = ctx->rate - ctx->pos;

        if (n > ctx->inlen) { n = ctx->inlen; }

        for (size_t i = 0; i < n; i++) { s[ctx->pos + i] ^= ctx->in[i]; }

        ctx->pos += n;
        ctx->in += n;
        ctx->inlen -= n;

        if (ctx->pos == ctx->rate)
        {
            keccak_f1600(ctx->state);
            ctx->pos = 0;
        }
    }
}

static void shake_finalize(struct shake_ctx *ctx)
{
    uint8_t *s = (uint8_t *)ctx->state;

    s[ctx->pos] ^= 0x1F;
    s[ctx->rate - 1] ^= 0x80;

    keccak_f1600(ctx->state);
}

static void shake_squeeze(struct shake_ctx *ctx)
{
    uint8_t *s = (uint8_t *)ctx->state;

    while (ctx->outlen > 0)
    {
        if (ctx->pos == ctx->rate)
        {
            keccak_f1600(ctx->state);
            ctx->pos = 0;
        }

        size_t n = ctx->rate - ctx->pos;

        if (n > ctx->outlen) { n = ctx->outlen; }

        for (size_t i = 0; i < n; i++) { ctx->out[i] = s[ctx->pos + i]; }

        ctx->pos += n;
        ctx->out += n;
        ctx->outlen -= n;
    }
}



